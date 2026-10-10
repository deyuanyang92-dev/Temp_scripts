#!/usr/bin/env python3
"""mito_trna.py - multi-tool tRNA annotation of mitochondrial genomes, GeSeq style.

Several tRNA predictors are run on each record and ALL their predictions are written into one
GenBank file for review (e.g. in Geneious), as GeSeq does with its third-party tRNA tools
(Tillich et al. 2017, Nucleic Acids Res 45:W6, doi:10.1093/nar/gkx391). A second GenBank keeps one
tRNA per locus for submission. The number of tRNAs is never imposed: it varies widely in animal
mtDNA (e.g. 1-2 in most anthozoans, 23-24 with duplications in some annelids), so every
prediction is judged on its own evidence.

Sources (any subset, --tools):
  mitos2    MITOS2 runmitos, tRNA module only (metazoan reference data, e.g. refseq89m)
  arwen     ARWEN (metazoan mitochondrial tRNA; recommended by GeSeq for metazoan mtDNA)
  aragorn   ARAGORN (-mt)
  trnascan  tRNAscan-SE 2 (-O organellar mode)
  input     tRNA features already present in the input GenBank (e.g. MitoZ)

Support is counted in independent method families (MITOS2; ARWEN+ARAGORN as one, same author and
algorithm; tRNAscan-SE; input annotation), not in tools.

Locus groups: predictions on the same strand overlapping >= 50 % of the shorter one, with the
same amino acid AND anticodon (so trnL1/trnL2 and trnS1/trnS2 stay separate). Overlapping
predictions with different anticodons form separate groups flagged as a conflict.

Selection (--select mitos2_or_support, default): a group is selected when MITOS2 found it, or when
>= --min_support (2) independent method families agree. Overlapping groups with different anticodons
(conflict): the one with stronger evidence (MITOS2, then more families) is kept, ties drop both.
--select mitos2: the input tRNAs (e.g. MitoZ) are replaced by the MITOS2 tRNAs (tools mitos2,input;
the input tRNAs stay visible in the review GenBank). Boundaries otherwise come from the first source
in --priority (MITOS2 first: most exact boundaries vs NCBI).

MITOS2 runs in fast mode by default (scripts/runmitos_fast.py): runmitos re-searches every tRNA
model it did not find with a sensitive cmsearch (--max, E 0.1, single CPU): 640 s vs 8 s on one
anthozoan (Metridium), plus spurious tRNAs; fast mode skips that pass (--mitos_mode full keeps it).

Outputs per input file: review/<name>.gbf (all predictions; /note = source, score, anticodon,
group, support), selected/<name>.gbf (selected tRNAs only; input for mito_curate.py -i out/selected),
and trna_report.tsv (one row per group).

Usage:
  python3 mito_trna.py -i gb_or_fasta -o out_dir --genetic_code 5 \
      --mitos_refdir DIR --mitos_ref refseq89m --mitos_bin ENV/bin --bin ENV/bin
"""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
ONE2THREE = {v: k for k, v in AA3.items()}
GB_EXT = {".gb", ".gbk", ".gbf", ".genbank", ".gbff"}
TIMEOUT = 3600


@dataclass
class Pred:
    source: str
    start: int            # 1-based inclusive
    end: int
    strand: int
    aa: str               # one-letter
    anticodon: str        # 5'->3', lowercase dna (e.g. 'tag'); '' if unknown
    score: str = ""
    group: int = -1
    extra: Dict[str, str] = field(default_factory=dict)

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def run(cmd: List[str], cwd: Optional[Path] = None) -> str:
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=TIMEOUT, cwd=cwd).stdout
    except (subprocess.TimeoutExpired, OSError) as exc:
        print(f"[mito_trna] {cmd[0]} failed: {exc}", file=sys.stderr)
        return ""


# ── predictors ───────────────────────────────────────────────────────────────
LASLETT = re.compile(r"m?tRNA-([A-Z][a-z]{2}|\?\?\?)\S*\s+(c?)\[(\d+),(\d+)\]\s+(\S+)?\s*\(([acgtnACGTN]{3})\)")


def parse_laslett(text: str, source: str) -> List[Pred]:
    """ARWEN / ARAGORN text output ('mtRNA-Leu  c[100,165]  33  (tag)'). A tRNA across the origin of a
    circular sequence is reported as [start,end] with start > end; it is kept with end = end + L
    (extended coordinates, written as a join in GenBank)."""
    out = []
    for m in LASLETT.finditer(text):
        aa = AA3.get(m.group(1), "")
        if aa:
            s, e = int(m.group(3)), int(m.group(4))
            out.append(Pred(source, s, e, -1 if m.group(2) else 1, aa, m.group(6).lower(), m.group(5) or "",
                            extra={"wrap": "1"} if s > e else {}))
    return out


def unwrap(preds: List[Pred], L: int) -> List[Pred]:
    """Origin-spanning predictions (start > end) -> extended coordinates end + L."""
    for p in preds:
        if p.start > p.end:
            p.end += L
    return preds


def arwen(fa: Path, code: int, b: Dict[str, str], tmp: Path, a) -> List[Pred]:
    return parse_laslett(run([b["arwen"], f"-gc{code}", "-mtx", "-w", str(fa)]), "ARWEN")


def aragorn(fa: Path, code: int, b: Dict[str, str], tmp: Path, a) -> List[Pred]:
    return parse_laslett(run([b["aragorn"], "-mt", f"-gc{code}", "-w", str(fa)]), "ARAGORN")


def trnascan(fa: Path, code: int, b: Dict[str, str], tmp: Path, a) -> List[Pred]:
    res = tmp / "trnascan.out"
    res.unlink(missing_ok=True)
    run([b["tRNAscan-SE"], "-O", "-q", "--brief", "-o", str(res), str(fa)])
    out = []
    for line in res.read_text().splitlines() if res.exists() else []:
        c = [x.strip() for x in line.split("\t")]
        if len(c) < 9 or c[4][:3] not in AA3:
            continue
        s, e = int(c[2]), int(c[3])
        out.append(Pred("tRNAscan-SE", min(s, e), max(s, e), 1 if s < e else -1, AA3[c[4][:3]],
                        c[5].lower().replace("u", "t"), c[8]))
    return out


MITOS_NAME = re.compile(r"^trn([A-Z])\d?\(([acgtnu]{3})\)", re.I)


def mitos2(fa: Path, code: int, b: Dict[str, str], tmp: Path, a) -> List[Pred]:
    """MITOS2, tRNA prediction only (proteins, rRNA, introns and OL skipped for speed)."""
    od = tmp / "mitos"
    shutil.rmtree(od, ignore_errors=True)
    od.mkdir()
    env_bin = str(Path(b["runmitos"]).parent)
    exe = [b["runmitos"]] if a.mitos_mode == "full" else \
        [str(Path(env_bin) / "python"), str(Path(__file__).resolve().parent / "scripts" / "runmitos_fast.py")]
    cmd = ["env", f"PATH={env_bin}:/usr/bin:/bin", *exe, "-i", str(fa), "-c", str(code), "-o", str(od),
           "-r", a.mitos_ref, "--refdir", a.mitos_refdir, "--prot", "0", "--rrna", "0", "--intron", "0", "--oril", "0"]
    run(cmd)
    bed = od / "result.bed"
    if not bed.exists():                              # failed run, not "no tRNA found"
        raise RuntimeError("MITOS2 produced no result.bed")
    out = []
    for line in bed.read_text().splitlines():
        c = line.split("\t")
        if len(c) < 6:
            continue
        m = MITOS_NAME.match(c[3])
        if m:
            out.append(Pred("MITOS2", int(c[1]) + 1, int(c[2]), -1 if c[5].strip() == "-" else 1, m.group(1).upper(),
                            m.group(2).lower().replace("u", "t"), c[4], extra={"mitos_name": c[3]}))
    return out


TOOLS = {"mitos2": mitos2, "arwen": arwen, "aragorn": aragorn, "trnascan": trnascan}


def from_input(rec) -> List[Pred]:
    """tRNA features of the input annotation (amino acid from /product, anticodon from /note or name)."""
    out = []
    for f in rec.features:
        if f.type != "tRNA":
            continue
        txt = " ".join(f.qualifiers.get("product", []) + f.qualifiers.get("gene", []) + f.qualifiers.get("note", []))
        m = re.search(r"tRNA-([A-Z][a-z]{2})", txt)
        aa = AA3.get(m.group(1), "") if m else ""
        if not aa:
            g = re.match(r"trn([A-Z])", f.qualifiers.get("gene", [""])[0])
            aa = g.group(1) if g else ""
        ac = re.search(r"\(([acgtuACGTU]{3})\)", txt)
        if aa:
            out.append(Pred("input", int(f.location.start) + 1, int(f.location.end), f.location.strand or 1, aa,
                            ac.group(1).lower().replace("u", "t") if ac else ""))
    return out


# ── grouping and selection ──────────────────────────────────────────────────
CIRC_L = 0                                            # length of the current (circular) record


def overlap(p: "Pred", q: "Pred") -> int:
    """Overlap on a circular sequence: origin-spanning predictions use extended coordinates (end > L),
    so compare also with one of them shifted by -L."""
    best = min(p.end, q.end) - max(p.start, q.start) + 1
    if CIRC_L:
        for dp, dq in ((-CIRC_L, 0), (0, -CIRC_L)):
            best = max(best, min(p.end + dp, q.end + dq) - max(p.start + dp, q.start + dq) + 1)
    return best


def group(preds: List[Pred]) -> List[List[Pred]]:
    groups: List[List[Pred]] = []
    for p in sorted(preds, key=lambda x: (x.start, x.end)):
        for g in groups:
            q = g[0]
            ov = overlap(p, q)
            same = p.aa == q.aa and (not p.anticodon or not q.anticodon or p.anticodon == q.anticodon)
            if p.strand == q.strand and same and ov >= 0.5 * min(p.length, q.length):
                g.append(p)
                break
        else:
            groups.append([p])
    for i, g in enumerate(groups):
        for p in g:
            p.group = i + 1
    return groups


def conflicts(groups: List[List[Pred]]) -> Dict[int, List[int]]:
    """Groups overlapping another group on the same strand with a different aa/anticodon."""
    out: Dict[int, List[int]] = {}
    for i, g in enumerate(groups):
        for j, h in enumerate(groups):
            if i == j or g[0].strand != h[0].strand:
                continue
            if overlap(g[0], h[0]) >= 0.5 * min(g[0].length, h[0].length):
                out.setdefault(i + 1, []).append(j + 1)
    return out


# Independent method families: ARWEN and ARAGORN share author and algorithm (Laslett & Canback),
# so their agreement is one line of evidence, not two (Metridium: both report trnF/trnS2 where NCBI
# and the covariance-model methods have no tRNA).
FAMILY = {"MITOS2": "MITOS2 (covariance models, metazoan mt tRNAs)", "ARWEN": "ARWEN/ARAGORN (heuristic)",
          "ARAGORN": "ARWEN/ARAGORN (heuristic)", "tRNAscan-SE": "tRNAscan-SE (covariance models)",
          "input": "input annotation"}


def families(g: List[Pred]) -> List[str]:
    return sorted({FAMILY.get(p.source, p.source) for p in g})


def strength(g: List[Pred], select: str) -> tuple:
    """Evidence rank used to settle conflicts: MITOS2 found it (mitos2_or_support mode), then families."""
    return (select == "mitos2_or_support" and any(p.source == "MITOS2" for p in g), len(families(g)))


def anticodon_of(g: List[Pred]) -> str:
    c = Counter(p.anticodon for p in g if p.anticodon)
    return c.most_common(1)[0][0] if c else ""


def choose(g: List[Pred], priority: List[str]) -> Pred:
    order = {s: i for i, s in enumerate(priority)}
    return min(g, key=lambda p: order.get(p.source, 99))


def trna_name(aa: str, ac: str) -> str:
    """trnL1/trnL2 and trnS1/trnS2 by anticodon (invertebrate convention: L1 = CUN (tag), L2 = UUR (taa);
    S1 = AGN (tct/gct), S2 = UCN (tga))."""
    if aa == "L":
        return "trnL1" if ac.endswith("ag") else "trnL2" if ac.endswith("aa") else "trnL"
    if aa == "S":
        return "trnS1" if ac.endswith("ct") else "trnS2" if ac.endswith("ga") else "trnS"
    return f"trn{aa}"


def feature(p: Pred, ac: str, note: str, L: int) -> SeqFeature:
    name = trna_name(p.aa, ac)
    loc: Any
    if p.end > L:                                     # across the origin of a circular sequence
        parts = [FeatureLocation(p.start - 1, L, strand=p.strand), FeatureLocation(0, p.end - L, strand=p.strand)]
        loc = CompoundLocation(parts if p.strand == 1 else parts[::-1])
    else:
        loc = FeatureLocation(p.start - 1, p.end, strand=p.strand)
    return SeqFeature(loc, type="tRNA",
                      qualifiers={"gene": [name + (f"({ac.replace('t', 'u')})" if ac else "")],   # MitoZ style trnN(guu)
                                  "product": [f"tRNA-{ONE2THREE[p.aa]}"], "note": [note]})


# ── driver ───────────────────────────────────────────────────────────────────
def find(name: str, *dirs: Optional[str]) -> Optional[str]:
    for d in dirs:
        if d and (Path(d) / name).exists():
            return str(Path(d) / name)
    return shutil.which(name)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="GenBank or FASTA file / directory")
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--genetic_code", type=int, default=5)
    p.add_argument("--tools", default=None,
                   help="default: mitos2,arwen,aragorn,trnascan,input; with --select mitos2: mitos2,input")
    p.add_argument("--priority", default="MITOS2,ARWEN,ARAGORN,tRNAscan-SE,input",
                   help="source whose boundaries are used for the selected tRNA of a group")
    p.add_argument("--select", choices=["mitos2_or_support", "support", "mitos2"], default="mitos2_or_support",
                   help="selected tRNAs: mitos2 = replace the input (e.g. MitoZ) tRNAs with the MITOS2 tRNAs "
                        "(MITOS2 boundaries and anticodons; a record whose MITOS2 run fails is kept unchanged); "
                        "mitos2_or_support (default) = found by MITOS2, or supported by >= "
                        "--min_support independent method families; support = families only. Benchmark "
                        "(39 annelid + 71 anthozoan/sponge/placozoan NCBI mitogenomes): MITOS2 alone had the "
                        "fewest false positives; ARWEN/ARAGORN alone gave hundreds in tRNA-poor genomes")
    p.add_argument("--min_support", type=int, default=2,
                   help="independent method families needed for a group to be selected (ARWEN+ARAGORN count once)")
    p.add_argument("--bin", default=None, help="directory with arwen, aragorn, tRNAscan-SE")
    p.add_argument("--mitos_bin", default=None, help="directory with runmitos (MITOS2 environment bin)")
    p.add_argument("--mitos_refdir", default=None, help="MITOS2 reference data directory (contains refseq89m)")
    p.add_argument("--mitos_ref", default="refseq89m")
    p.add_argument("--mitos_mode", choices=["fast", "full"], default="fast",
                   help="fast (default): MITOS2 without its sensitive-mode fallback (scripts/runmitos_fast.py); "
                        "full: unmodified runmitos (much slower when tRNAs are missing, e.g. anthozoans)")
    a = p.parse_args(argv)

    if a.tools is None:
        a.tools = "mitos2,input" if a.select == "mitos2" else "mitos2,arwen,aragorn,trnascan,input"
    want = [t.strip() for t in a.tools.split(",") if t.strip()]
    found = {n: find(n, a.bin, a.mitos_bin) for n in ("arwen", "aragorn", "tRNAscan-SE", "runmitos")}
    need = {"mitos2": "runmitos", "arwen": "arwen", "aragorn": "aragorn", "trnascan": "tRNAscan-SE"}
    for t in list(want):
        if t in need and not found[need[t]]:
            print(f"[mito_trna] {t}: {need[t]} not found - skipped", file=sys.stderr)
            want.remove(t)
        if t == "mitos2" and t in want and not a.mitos_refdir:
            print("[mito_trna] mitos2: --mitos_refdir not given - skipped", file=sys.stderr)
            want.remove(t)
    if a.select == "mitos2" and "mitos2" not in want:
        print("[mito_trna] --select mitos2 needs MITOS2 (runmitos via --mitos_bin, and --mitos_refdir)", file=sys.stderr)
        return 2
    b: Dict[str, str] = {k: v for k, v in found.items() if v}
    src, out = Path(a.input), Path(a.out)
    files = sorted(x for x in src.rglob("*") if x.is_file() and x.suffix.lower() in GB_EXT | {".fa", ".fasta", ".fna"}) \
        if src.is_dir() else [src]
    if a.select == "mitos2":
        # MITOS2 runs on the sequence of each input GenBank record, so its coordinates always match that
        # record; the input must therefore be the final, reoriented annotation (batch_mitoz 04.final_gb),
        # never a FASTA or the round-1 MitoZ output from before reorientation
        bad = [x for x in files if x.suffix.lower() not in GB_EXT or any(
            part.startswith(("00.", "01.mitoz_anno1")) for part in x.parts)]
        if bad:
            print(f"[mito_trna] --select mitos2 needs the final reoriented GenBank (e.g. 04.final_gb/), "
                  f"not: {', '.join(map(str, bad[:3]))}", file=sys.stderr)
            return 2
    out.mkdir(parents=True, exist_ok=True)
    priority = [x.strip() for x in a.priority.split(",")]
    report: List[Dict[str, str]] = []
    for f in files:
        fmt = "genbank" if f.suffix.lower() in GB_EXT else "fasta"
        recs = list(SeqIO.parse(str(f), fmt))
        review, selected = [], []
        for rec in recs:
            rec.annotations.setdefault("molecule_type", "DNA")
            with tempfile.TemporaryDirectory() as td:
                tmp = Path(td)
                fa = tmp / "g.fa"
                fa.write_text(f">{rec.id}\n{rec.seq}\n")
                preds: List[Pred] = []
                failed = []
                for t in want:
                    if t == "input":
                        preds += from_input(rec)
                        continue
                    try:
                        preds += TOOLS[t](fa, a.genetic_code, b, tmp, a)
                    except RuntimeError as exc:
                        failed.append(t)
                        print(f"[mito_trna] {rec.id}: {t} failed ({exc})", file=sys.stderr)
            if a.select == "mitos2" and "mitos2" in failed:
                # never drop the existing tRNAs because MITOS2 crashed: record written unchanged
                print(f"[mito_trna] {rec.id}: MITOS2 failed - input tRNAs kept unchanged", file=sys.stderr)
                review.append(rec)
                selected.append(rec)
                report.append({"file": f.name, "record": rec.id, "group": "", "tRNA": "", "anticodon": "", "strand": "",
                               "start": "", "end": "", "boundary_source": "", "n_families": "", "sources": "",
                               "coords_by_source": "", "conflict": "", "selected": "MITOS2 failed; input kept"})
                continue
            L = len(rec.seq)
            global CIRC_L
            CIRC_L = L if rec.annotations.get("topology", "circular") == "circular" else 0
            preds = unwrap(preds, L)
            groups = group(preds)
            conf = conflicts(groups)
            sources_run = sorted({x for x in ("MITOS2" if "mitos2" in want else "", "ARWEN" if "arwen" in want else "",
                                              "ARAGORN" if "aragorn" in want else "",
                                              "tRNAscan-SE" if "trnascan" in want else "",
                                              "input" if "input" in want else "") if x})
            # old tRNAs and their gene features are replaced (MitoZ writes a gene feature per tRNA)
            trna_genes = {g for x in rec.features if x.type == "tRNA" for g in x.qualifiers.get("gene", [])}
            others = [x for x in rec.features if x.type != "tRNA" and not (
                x.type == "gene" and (set(x.qualifiers.get("gene", [])) & trna_genes
                                      or x.qualifiers.get("gene", [""])[0].lower().startswith("trn")))]
            rev_feats, sel_feats = list(others), list(others)
            for gi, g in enumerate(groups, 1):
                ac = anticodon_of(g)
                srcs = sorted({x.source for x in g})
                fams = families(g)
                n_fam_run = len({FAMILY.get(x, x) for x in sources_run})
                support = (f"{len(fams)}/{n_fam_run} independent method families; sources: {', '.join(srcs)}")
                flag = f"; CONFLICT with group(s) {conf[gi]}" if gi in conf else ""
                for x in g:
                    rev_feats.append(feature(x, ac, f"source: {x.source}; score {x.score or 'NA'}; anticodon "
                                                    f"{x.anticodon or 'NA'}; group {gi}; support {support}{flag}", L))
                best = choose(g, priority)
                # a conflict removes a group only when a conflicting group has equal or stronger
                # evidence (ARWEN-only hits overlapping a MITOS2 tRNA must not remove it)
                beaten = any(strength(groups[j - 1], a.select) >= strength(g, a.select) for j in conf.get(gi, []))
                if a.select == "mitos2":
                    keep = "MITOS2" in srcs
                    best = next(x for x in g if x.source == "MITOS2") if keep else best
                elif a.select == "mitos2_or_support":
                    keep = ("MITOS2" in srcs or len(fams) >= a.min_support) and not beaten
                else:
                    keep = len(fams) >= a.min_support and not beaten
                if keep:
                    ac_sel = (best.anticodon or ac) if a.select == "mitos2" else ac
                    tf = feature(best, ac_sel, f"tRNA by mito_trna.py: boundaries from {best.source}; "
                                               f"support {support}", L)
                    sel_feats += [SeqFeature(tf.location, type="gene", qualifiers={"gene": tf.qualifiers["gene"]}), tf]
                report.append({"file": f.name, "record": rec.id, "group": str(gi), "tRNA": trna_name(best.aa, ac),
                               "anticodon": ac, "strand": "+" if best.strand == 1 else "-",
                               "start": str(best.start), "end": str(best.end), "boundary_source": best.source,
                               "n_families": str(len(fams)), "sources": ",".join(srcs),
                               "coords_by_source": "; ".join(f"{x.source}:{x.start}..{x.end}" for x in g),
                               "conflict": ",".join(map(str, conf.get(gi, []))), "selected": "yes" if keep else "no"})
            key = lambda x: (int(x.location.start), x.type != "source", x.type != "gene")  # noqa: E731
            r1, r2 = rec[:], rec[:]
            r1.annotations, r2.annotations = dict(rec.annotations), dict(rec.annotations)
            r1.features, r2.features = sorted(rev_feats, key=key), sorted(sel_feats, key=key)
            review.append(r1)
            selected.append(r2)
        for sub, recs_out in (("review", review), ("selected", selected)):
            (out / sub).mkdir(exist_ok=True)
            SeqIO.write(recs_out, str(out / sub / f"{f.stem}.gbf"), "genbank")
    cols = ["file", "record", "group", "tRNA", "anticodon", "strand", "start", "end", "boundary_source",
            "n_families", "sources", "coords_by_source", "conflict", "selected"]
    with open(out / "trna_report.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(report)
    sel = sum(r["selected"] == "yes" for r in report)
    rule = "MITOS2 tRNAs replace input tRNAs" if a.select == "mitos2" else \
        ("found by MITOS2 or " if a.select == "mitos2_or_support" else "") + \
        f">= {a.min_support} independent method families"
    print(f"[mito_trna] {len(files)} file(s); tools {', '.join(want)}; {len(report)} tRNA loci, {sel} selected "
          f"({rule}{'' if a.select == 'mitos2' else '; conflicts -> stronger evidence'}) -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
