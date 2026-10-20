"""Evaluate benchmark runs against the NCBI (NC_) annotation.

Features: the 13 PCGs, rrnS/rrnL, tRNAs, introns (others, e.g. ORFs and control region, are ignored).
A prediction matches an NCBI feature of the same kind when it overlaps >= 50 % of it; among the
overlapping ones the same name + strand is preferred. Classes per NCBI feature: exact, le3 (both ends
within 3 nt), gt3, strand_error, name_error, name_unverifiable (NCBI gives no anticodon for trnL/trnS),
missed; unmatched predictions: fp. Obvious errors (spec §5.2): missing PCG/rRNA, internal stop,
strand error, duplicate gene, name error, boundary > 30 nt, table2asn ERROR/REJECT, failed/timeout.

Outputs in <root>/eval/<split>/: features.tsv, pcg_checks.tsv, errors.tsv, summary.tsv, timing.tsv.

Usage: python3 mito_bench_eval.py --set sets/mito_bench_2026-10.tsv --root projects/mito_benchmark
         [--split dev] [--tools ours,mitoz,mitos2] [--refseq_cutoff YYYY-MM-DD] [--table2asn]
"""
import argparse
import csv
import json
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, NamedTuple

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from mito_names import PCG, gene_name, trna_name  # noqa: E402

MIN_OVERLAP = 0.5
NEAR = 3
OBVIOUS_NT = 30
AMBIGUOUS_AA = {"L", "S"}              # two tRNA genes each; told apart only by anticodon
FOUND = ("exact", "le3", "gt3", "name_unverifiable")


class Feat(NamedTuple):
    kind: str          # PCG, rRNA, tRNA, intron
    name: str          # cox1 / rrnL / amino acid letter / intron
    strand: int
    start: int         # 0-based, end exclusive (span of all parts)
    end: int
    anticodon: str


def norm_features(rec: SeqRecord) -> List[Feat]:
    out = []
    for f in rec.features:
        strand = f.location.strand or 1
        s, e = int(f.location.start), int(f.location.end)
        if f.type in ("CDS", "rRNA"):
            g = gene_name(f.qualifiers, f.type)
            if g:
                out.append(Feat("PCG" if g in PCG else "rRNA", g, strand, s, e, ""))
        elif f.type == "tRNA":
            aa, ac = trna_name(f.qualifiers)
            if aa:
                out.append(Feat("tRNA", aa, strand, s, e, ac))
        elif f.type == "intron":
            out.append(Feat("intron", "intron", strand, s, e, ""))
    return out


def _overlap(a: Feat, b: Feat) -> int:
    return max(0, min(a.end, b.end) - max(a.start, b.start))


def _same_name(t: Feat, p: Feat) -> str:
    """'yes', 'no' or 'unverifiable' (tRNA of an ambiguous amino acid without anticodon on one side)."""
    if t.name != p.name:
        return "no"
    if t.kind != "tRNA" or t.name not in AMBIGUOUS_AA:
        return "yes"
    if not t.anticodon or not p.anticodon:
        return "unverifiable"
    # anticodon bases 2-3 pair with codon positions 2-1: they fix the codon box (Ser AGN vs UCN, Leu CUN vs UUR);
    # the wobble base 1 may differ between annotations of the same gene (e.g. NCBI 'act' vs 'gct')
    return "yes" if t.anticodon[1:] == p.anticodon[1:] else "no"


def match(truth: List[Feat], pred: List[Feat]) -> List[dict]:
    rows, used = [], set()
    for t in sorted(truth, key=lambda x: (x.kind, x.start)):
        cands = [(i, p) for i, p in enumerate(pred) if i not in used and p.kind == t.kind
                 and _overlap(t, p) >= MIN_OVERLAP * (t.end - t.start)]
        row: dict = {"kind": t.kind, "name": t.name, "truth": f"{t.start + 1}..{t.end}{'+' if t.strand > 0 else '-'}",
               "pred": "", "delta": "", "class": "missed"}
        if cands:
            def rank(c):
                same = _same_name(t, c[1]) != "no"
                return (not (same and c[1].strand == t.strand), not same, c[1].strand != t.strand, -_overlap(t, c[1]))
            i, p = min(cands, key=rank)
            used.add(i)
            same = _same_name(t, p)
            delta = max(abs(p.start - t.start), abs(p.end - t.end))
            row.update(pred=f"{p.start + 1}..{p.end}{'+' if p.strand > 0 else '-'}", delta=delta)
            if same == "no":
                row["class"] = "name_error"
            elif p.strand != t.strand:
                row["class"] = "strand_error"
            elif same == "unverifiable":
                row["class"] = "name_unverifiable"
            else:
                row["class"] = "exact" if delta == 0 else "le3" if delta <= NEAR else "gt3"
        rows.append(row)
    for i, p in enumerate(pred):
        if i not in used:
            rows.append({"kind": p.kind, "name": p.name, "truth": "",
                         "pred": f"{p.start + 1}..{p.end}{'+' if p.strand > 0 else '-'}", "delta": "", "class": "fp"})
    return rows


def duplicates(pred: List[Feat]) -> List[str]:
    """PCG / rRNA names annotated more than once at non-overlapping places."""
    by: Dict[str, List[Feat]] = defaultdict(list)
    for p in pred:
        if p.kind in ("PCG", "rRNA"):
            by[p.name].append(p)
    return sorted(n for n, fs in by.items() if any(_overlap(a, b) == 0 for i, a in enumerate(fs) for b in fs[i + 1:]))


def pcg_checks(rec: SeqRecord, code: int) -> List[dict]:
    """internal_stop / bad_start / bad_stop per standard PCG; incomplete stops T / TA are accepted."""
    table = CodonTable.unambiguous_dna_by_id[code]
    out = []
    for f in rec.features:
        if f.type != "CDS" or not gene_name(f.qualifiers, "CDS"):
            continue
        if "<" in str(f.location) or ">" in str(f.location):
            continue
        g = gene_name(f.qualifiers, "CDS")
        seq = str(f.extract(rec.seq)).upper()[int(f.qualifiers.get("codon_start", ["1"])[0]) - 1:]
        rest = len(seq) % 3
        body = seq[:len(seq) - rest]
        codons = [body[i:i + 3] for i in range(0, len(body), 3)]
        if not codons:
            continue
        stop_ok = rest == 0 and codons[-1] in table.stop_codons
        if rest:
            stop_ok = seq[-rest:] in ("T", "TA") or bool(f.qualifiers.get("transl_except"))
        inner = codons[:-1] if rest == 0 else codons
        if any(c in table.stop_codons for c in inner[1:]):
            out.append({"gene": g, "check": "internal_stop"})
        if codons[0] not in table.start_codons:
            out.append({"gene": g, "check": "bad_start"})
        if not stop_ok:
            out.append({"gene": g, "check": "bad_stop"})
    return out


def obvious_errors(rows: List[dict], pcg_rows: List[dict], status: dict) -> List[dict]:
    st = status.get("status")
    if st == "unsupported":
        return []
    if st != "ok":
        return [{"error": st, "gene": "", "detail": status.get("error", "")[:200]}]
    out = []
    for r in rows:
        c = r["class"]
        if c == "missed" and r["kind"] in ("PCG", "rRNA"):
            out.append({"error": "missing_gene", "gene": r["name"], "detail": r.get("truth", "")})
        elif c in ("strand_error", "name_error"):
            out.append({"error": c, "gene": r["name"], "detail": f"{r.get('truth', '')} vs {r.get('pred', '')}"})
        elif c == "gt3" and int(r["delta"]) > OBVIOUS_NT:
            out.append({"error": "boundary_gt30", "gene": r["name"], "detail": f"{r.get('truth', '')} vs {r.get('pred', '')}"})
    out += [{"error": "internal_stop", "gene": p["gene"], "detail": ""} for p in pcg_rows if p["check"] == "internal_stop"]
    return out


# ── driver ───────────────────────────────────────────────────────────────────
def table2asn_errors(files: Dict[str, Path], work: Path) -> Dict[str, int]:
    """accession -> number of table2asn ERROR/REJECT messages (via ncbi_validate.py)."""
    src = work / "in"
    src.mkdir(parents=True, exist_ok=True)
    for acc, p in files.items():
        (src / f"{acc}.gbf").write_text(p.read_text())
    subprocess.run([sys.executable, str(HERE.parent / "ncbi_validate.py"), "-i", str(src), "-o", str(work / "out")],
                   capture_output=True, text=True)
    n: Counter = Counter()
    val = work / "out" / "validation.tsv"
    if val.exists():
        for r in csv.DictReader(val.open(encoding="utf-8"), delimiter="\t"):
            if r["severity"] in ("ERROR", "REJECT"):
                n[r["record"].split(".")[0]] += 1
    return {acc: n.get(acc.split(".")[0], 0) for acc in files}


def strata(row: dict, cutoff: str) -> List[tuple]:
    out = [("all", "all"), ("phylum", row["phylum"]), ("evidence_tier", row["evidence_tier"]),
           ("annotation_method", row["annotation_method"]),
           ("structure_stress", "yes" if row["structure_flags"] else "no")]
    if cutoff:
        out.append(("post_cutoff", "yes" if row["refseq_created"] > cutoff else "no"))
    return out


def write_tsv(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, keys, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--set", required=True)
    p.add_argument("--root", required=True)
    p.add_argument("--split", default="dev")
    p.add_argument("--tools", default="ours,mitoz,mitos2")
    p.add_argument("--refseq_cutoff", default="", help="records created after this date form the post_cutoff stratum")
    p.add_argument("--table2asn", action="store_true", help="also count NCBI validator ERROR/REJECT per record")
    a = p.parse_args(argv)
    root = Path(a.root)
    out = root / "eval" / a.split
    out.mkdir(parents=True, exist_ok=True)
    rows = [r for r in csv.DictReader(open(a.set, encoding="utf-8"), delimiter="\t")
            if a.split == "all" or r["split"] == a.split]
    tools = a.tools.split(",")
    feats, checks, errors, timing = [], [], [], []
    status: Dict[tuple, dict] = {}
    finals: Dict[str, Dict[str, Path]] = defaultdict(dict)
    for r in rows:
        acc, code = r["accession"], int(r["genetic_code"])
        ncbi = SeqIO.read(str(root / "gb" / f"{acc}.gb"), "genbank")
        truth = norm_features(ncbi)
        checks += [{"tool": "ncbi", "accession": acc, **c} for c in pcg_checks(ncbi, code)]
        for t in tools:
            sp = root / "runs" / t / acc / "status.json"
            st = json.loads(sp.read_text()) if sp.exists() else {"status": "not_run", "seconds": 0, "error": ""}
            status[(t, acc)] = st
            timing.append({"tool": t, "accession": acc, "status": st["status"], "seconds": st.get("seconds", 0),
                           "steps": json.dumps(st.get("step_seconds", {}))})
            if st["status"] == "not_run":
                continue
            m, pc, dup = [], [], []
            if st["status"] == "ok":
                fp = root / "runs" / t / acc / "final.gbf"
                finals[t][acc] = fp
                pred_rec = SeqIO.read(str(fp), "genbank")
                pred = norm_features(pred_rec)
                m = match(truth, pred)
                pc = pcg_checks(pred_rec, code)
                dup = duplicates(pred)
            feats += [{"tool": t, "accession": acc, **x} for x in m]
            checks += [{"tool": t, "accession": acc, **c} for c in pc]
            errs = obvious_errors(m, pc, st) + [{"error": "duplicate_gene", "gene": g, "detail": ""} for g in dup]
            errors += [{"tool": t, "accession": acc, "phylum": r["phylum"], **e} for e in errs]
    if a.table2asn:
        for t in tools:
            for acc, n in table2asn_errors(finals[t], out / "table2asn" / t).items():
                if n:
                    errors.append({"tool": t, "accession": acc, "phylum": next(x["phylum"] for x in rows if x["accession"] == acc),
                                   "error": "table2asn_error", "gene": "", "detail": f"{n} ERROR/REJECT"})
    # summary: per tool x kind x stratum; 'common_ok' = records every tool finished
    common = {r["accession"] for r in rows if all(status.get((t, r["accession"]), {}).get("status") == "ok" for t in tools)}
    by_acc = {r["accession"]: r for r in rows}
    err_acc = {(e["tool"], e["accession"]) for e in errors}
    summ: Dict[tuple, Counter] = defaultdict(Counter)
    for x in feats:
        r = by_acc[x["accession"]]
        st_list = strata(r, a.refseq_cutoff) + ([("common_ok", "yes")] if x["accession"] in common else [])
        for s in st_list:
            c = summ[(x["tool"], x["kind"]) + s]
            c[x["class"]] += 1
            if x["class"] != "fp":
                c["n_truth"] += 1
    for (t, acc), st in status.items():
        r = by_acc[acc]
        st_list = strata(r, a.refseq_cutoff) + ([("common_ok", "yes")] if acc in common else [])
        for s in st_list:
            c = summ[(t, "record") + s]
            c["records"] += 1
            c[st["status"]] += 1
            if st["status"] == "ok" and (t, acc) in err_acc:
                c["with_obvious_error"] += 1
    srows = []
    for (t, kind, sname, sval), c in sorted(summ.items()):
        row = {"tool": t, "kind": kind, "stratum": sname, "value": sval}
        if kind == "record":
            row.update({k: c[k] for k in ("records", "ok", "failed", "timeout", "unsupported", "with_obvious_error")})
        else:
            found = sum(c[k] for k in FOUND)
            row.update({"n_truth": c["n_truth"], "found": found,
                        "recall": round(found / c["n_truth"], 4) if c["n_truth"] else "",
                        **{k: c[k] for k in ("exact", "le3", "gt3", "strand_error", "name_error",
                                             "name_unverifiable", "missed", "fp")}})
        srows.append(row)
    for name, data in (("features", feats), ("pcg_checks", checks), ("errors", errors), ("summary", srows),
                       ("timing", timing)):
        write_tsv(out / f"{name}.tsv", data)
    print(f"[eval] {len(rows)} records, {len(tools)} tools -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
