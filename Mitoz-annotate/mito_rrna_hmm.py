#!/usr/bin/env python3
"""mito_rrna_hmm.py - last-resort 12S/16S search with per-phylum profile HMMs (MitoM).

Most genomes are annotated well by the cores (MitoZ / MITOS2 covariance models); this step runs only for
an rRNA that is still missing after them and the MITOS2 rescue. It searches the genome with an nhmmer
profile HMM built from the NCBI rRNAs of the record's own phylum, as EMMA does with its (vertebrate)
HMMs (knowledge/tools/emma.md). 2026-10-21 dev benchmark, search level, leave-genus-out
(knowledge/topics/rrna_annotation_methods.md section 6): MITOS2 CM 166 -> CM then per-phylum HMM 175 of 186
NCBI rRNAs; nine divergent phyla 30 -> 38 of 48 (both ctenophores, nematodes, placozoan, acoel).

Rules:
  - only rRNAs absent from the record are searched; the phylum is taken from the record's taxonomy
    (or --phylum); phyla with fewer than MIN_REFS reference rRNAs have no model and are skipped;
  - hits of one model are joined only when they are pieces of one gene: same strand, disjoint and
    consecutive model columns in the same order as in the genome, <= MAX_GAP nt apart. Hits over the
    same model columns are copies (mitochondrial genes can be duplicated) and are never joined;
  - the joined span is added when it is >= MIN_RRNA nt, <= MAX_LEN_FACTOR x the phylum median length and
    does not lie mostly (> 50 %) on an existing CDS / rRNA; ends are set later by mito_curate.py R1.

Usage:
  python3 mito_rrna_hmm.py build --gb refs/*.gb -o models/rrna_hmm [--exclude accessions.txt]
  python3 mito_rrna_hmm.py search -i 04.final_gb/ [--lib models/rrna_hmm] [--phylum Ctenophora]
build: one GenBank file per phylum (file name = phylum), e.g. from benchmark/fetch_refs-style NCBI queries.
search writes the GenBank files in place and rrna_hmm_report.tsv next to them.

HMMER: Eddy 2011, PLoS Comput Biol 7:e1002195 (doi:10.1371/journal.pcbi.1002195); nhmmer: Wheeler &
Eddy 2013, Bioinformatics 29:2487 (doi:10.1093/bioinformatics/btt403); MAFFT: Katoh & Standley 2013,
Mol Biol Evol 30:772 (doi:10.1093/molbev/mst010).
"""
import argparse
import csv
import random
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Callable, Dict, Iterable, List, NamedTuple, Optional, Set

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402

LIB = Path(__file__).resolve().parent / "models" / "rrna_hmm"
GENES = {"rrnS": ("s-rRNA", "12S ribosomal RNA"), "rrnL": ("l-rRNA", "16S ribosomal RNA")}
NAMES = {"rrnS": {"rrns", "s-rrna", "12s", "rns", "12s rrna"}, "rrnL": {"rrnl", "l-rrna", "16s", "rnl", "16s rrna"}}
MIN_REFS, MAX_ALN, MIN_RRNA, MAX_GAP, MAX_LEN_FACTOR, EVALUE = 3, 60, 300, 300, 1.5, 1e-3
NOTE = ("not found by the annotation cores; added from a per-phylum profile HMM ({}, {} reference rRNAs) "
        "by mito_rrna_hmm.py - check boundaries manually")


class Hit(NamedTuple):
    start: int      # 0-based genome start (envelope)
    end: int        # genome end, exclusive
    strand: int
    hmm_from: int   # model columns, 1-based
    hmm_to: int
    score: float


def gene_of(f: SeqFeature) -> Optional[str]:
    """rrnS / rrnL for an rRNA feature named in MitoZ, NCBI or free style; None otherwise."""
    s = " ".join(f.qualifiers.get("gene", []) + f.qualifiers.get("product", [])).lower()
    for g, keys in NAMES.items():
        if any(k in s for k in keys) or ("small" in s if g == "rrnS" else "large" in s):
            return g
    return None


def missing(rec: SeqRecord) -> List[str]:
    have = {gene_of(f) for f in rec.features if f.type == "rRNA"}
    return [g for g in GENES if g not in have]


def phylum_of(rec: SeqRecord, lib_phyla: Iterable[str], override: Optional[str] = None) -> Optional[str]:
    lib = {p.lower(): p for p in lib_phyla}
    if override:
        return lib.get(override.lower())
    tax = rec.annotations.get("taxonomy")
    taxonomy: List[str] = [str(t) for t in tax] if isinstance(tax, list) else []
    return next((lib[t.lower()] for t in taxonomy if t.lower() in lib), None)


def _joinable(span: Hit, h: Hit) -> bool:
    if h.strand != span.strand:
        return False
    after = h.start >= span.end - 10                         # h lies 3' of span in the genome (plus strand)
    before = h.end <= span.start + 10
    if not (after or before):
        return False
    gap = h.start - span.end if after else span.start - h.end
    if gap > MAX_GAP:
        return False
    later_in_model = h.hmm_from >= span.hmm_to - 10
    earlier_in_model = h.hmm_to <= span.hmm_from + 10
    if span.strand == 1:
        return (after and later_in_model) or (before and earlier_in_model)
    return (after and earlier_in_model) or (before and later_in_model)


def best_span(hits: List[Hit]) -> Optional[Hit]:
    """Best hit, extended by the hits that are other pieces of the same gene (see module rules)."""
    if not hits:
        return None
    rest = sorted(hits, key=lambda h: -h.score)
    span = rest.pop(0)
    grown = True
    while grown:
        grown = False
        for h in rest:
            if _joinable(span, h):
                span = Hit(min(span.start, h.start), max(span.end, h.end), span.strand,
                           min(span.hmm_from, h.hmm_from), max(span.hmm_to, h.hmm_to), span.score + h.score)
                rest.remove(h)
                grown = True
                break
    return span


def read_manifest(lib: Path) -> Dict[tuple, dict]:
    f = lib / "manifest.tsv"
    if not f.exists():
        return {}
    return {(r["phylum"], r["gene"]): r for r in csv.DictReader(open(f), delimiter="\t")}


def nhmmer(hmm: Path, fa: Path) -> List[Hit]:
    exe = mito_env.find("nhmmer")
    if not exe:
        raise RuntimeError("nhmmer not found (HMMER 3; conda env mitos)")
    with tempfile.TemporaryDirectory() as td:
        tbl = Path(td) / "hits.tbl"
        subprocess.run([exe, "--dna", "--cpu", "2", "-E", str(EVALUE), "--tblout", str(tbl), "-o", "/dev/null",
                        str(hmm), str(fa)], check=True, capture_output=True)
        out = []
        for line in tbl.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            c = line.split()
            a, b, st = int(c[8]), int(c[9]), c[11]             # envelope from / to, strand
            s, e = (a - 1, b) if st == "+" else (b - 1, a)
            out.append(Hit(s, e, 1 if st == "+" else -1, int(c[4]), int(c[5]), float(c[13])))
        return out


def _positions(f: SeqFeature, n: int) -> Set[int]:
    return {i % n for p in f.location.parts for i in range(int(p.start), int(p.end))}


def search_record(rec: SeqRecord, lib: Path = LIB, phylum: Optional[str] = None,
                  run: Optional[Callable[[Path, Path], List[Hit]]] = None) -> List[dict]:
    """Search the rRNAs missing from rec; add the accepted ones in place. One report row per searched gene."""
    want = missing(rec)
    if not want:
        return []
    man = read_manifest(lib)
    phy = phylum_of(rec, {p for p, _ in man}, phylum)
    row = {"record": rec.id, "phylum": phy or ""}
    if not phy:
        return [dict(row, gene=GENES[g][0], location="", status="skipped: no library for this phylum") for g in want]
    run = run or nhmmer
    n = len(rec)
    taken = [_positions(f, n) for f in rec.features if f.type in ("CDS", "rRNA")]
    rows = []
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "q.fa"
        fa.write_text(f">{rec.id}\n{rec.seq}\n")
        for g in want:
            name, product = GENES[g]
            m = man.get((phy, g))
            if not m:
                rows.append(dict(row, gene=name, location="", status="skipped: no model for this gene"))
                continue
            span = best_span(run(lib / f"{phy}.{g}.hmm", fa))
            if span is None:
                rows.append(dict(row, gene=name, location="", status="not found"))
                continue
            loc = f"{span.start + 1}..{span.end}{'+' if span.strand == 1 else '-'}"
            length, pos = span.end - span.start, set(range(span.start, span.end))
            status = "added"
            if length < MIN_RRNA:
                status = f"rejected: short ({length} nt < {MIN_RRNA})"
            elif length > MAX_LEN_FACTOR * float(m["median_len"]):
                status = f"rejected: long ({length} nt > {MAX_LEN_FACTOR} x median {m['median_len']})"
            elif any(len(pos & t) > 0.5 * length for t in taken):
                status = "rejected: lies on an existing gene"
            rows.append(dict(row, gene=name, location=loc, status=status))
            if status != "added":
                continue
            fl = FeatureLocation(span.start, span.end, span.strand)
            q = {"gene": [name], "product": [product], "note": [NOTE.format(phy, m["n_used"])]}
            rec.features += [SeqFeature(fl, type="gene", qualifiers={"gene": [name]}), SeqFeature(fl, type="rRNA", qualifiers=q)]
            taken.append(pos)
    if any(r["status"] == "added" for r in rows):
        src = [f for f in rec.features if f.type == "source"]
        rec.features = src + sorted((f for f in rec.features if f.type != "source"),
                                    key=lambda f: (int(f.location.start), f.type != "gene"))
    return rows


def search_dir(final_dir: Path, lib: Path = LIB, phylum: Optional[str] = None,
               run: Optional[Callable[[Path, Path], List[Hit]]] = None) -> List[dict]:
    """All *.gbf / *.gb in final_dir, rewritten in place when an rRNA was added."""
    rows = []
    for f in sorted(list(final_dir.glob("*.gbf")) + list(final_dir.glob("*.gb"))):
        recs = list(SeqIO.parse(str(f), "genbank"))
        changed = False
        for r in recs:
            rs = search_record(r, lib, phylum, run)
            changed |= any(x["status"] == "added" for x in rs)
            rows += [dict(x, file=f.name) for x in rs]
        if changed:
            SeqIO.write(recs, str(f), "genbank")
    return rows


def _ref_rrnas(gb: Path, exclude: Set[str]) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {g: [] for g in GENES}
    for r in SeqIO.parse(str(gb), "genbank"):
        if r.id.split(".")[0] in exclude:
            continue
        for f in r.features:
            g = gene_of(f) if f.type == "rRNA" else None
            if not g:
                continue
            s = str(f.extract(r.seq)).upper()
            if 300 <= len(s) <= 2500 and s.count("N") <= 0.01 * len(s):
                out[g].append(s)
    return out


def build(gbs: List[Path], out: Path, exclude: Optional[Set[str]] = None, cpu: int = 4) -> List[dict]:
    """One HMM per phylum (GenBank file stem) and rRNA: MAFFT alignment of <= MAX_ALN random references,
    hmmbuild --dna. Writes out/<phylum>.<rrnS|rrnL>.hmm and out/manifest.tsv."""
    mafft, hmmbuild = mito_env.find("mafft"), mito_env.find("hmmbuild")
    if not (mafft and hmmbuild):
        raise RuntimeError("mafft and hmmbuild are needed to build the library")
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for gb in sorted(gbs):
        for g, seqs in _ref_rrnas(gb, exclude or set()).items():
            if len(seqs) < MIN_REFS:
                continue
            sub = seqs if len(seqs) <= MAX_ALN else random.Random(1).sample(seqs, MAX_ALN)
            with tempfile.TemporaryDirectory() as td:
                fa, aln = Path(td) / "in.fa", Path(td) / "aln.fa"
                fa.write_text("".join(f">r{i}\n{s}\n" for i, s in enumerate(sub)))
                aln.write_text(subprocess.run([mafft, "--auto", "--thread", str(cpu), "--quiet", str(fa)],
                                              capture_output=True, text=True, check=True).stdout)
                subprocess.run([hmmbuild, "--dna", "--cpu", str(cpu), "-n", f"{gb.stem}.{g}",
                                str(out / f"{gb.stem}.{g}.hmm"), str(aln)], capture_output=True, check=True)
            rows.append({"phylum": gb.stem, "gene": g, "n_ref": len(seqs), "n_used": len(sub),
                         "median_len": int(statistics.median(len(s) for s in seqs))})
    old = read_manifest(out)
    old.update({(r["phylum"], r["gene"]): r for r in rows})
    with open(out / "manifest.tsv", "w") as h:
        w = csv.DictWriter(h, fieldnames=["phylum", "gene", "n_ref", "n_used", "median_len"], delimiter="\t")
        w.writeheader()
        w.writerows(sorted(old.values(), key=lambda r: (r["phylum"], r["gene"])))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="build per-phylum rRNA HMMs from reference GenBank files")
    b.add_argument("--gb", nargs="+", type=Path, required=True, help="one GenBank file per phylum (stem = phylum)")
    b.add_argument("-o", "--out", type=Path, default=LIB)
    b.add_argument("--exclude", type=Path, help="file with accessions to leave out (one per line)")
    b.add_argument("--cpu", type=int, default=4)
    s = sub.add_parser("search", help="search missing rRNAs in final GenBank files (in place)")
    s.add_argument("-i", "--input", type=Path, required=True, help="directory of *.gbf / *.gb")
    s.add_argument("--lib", type=Path, default=LIB)
    s.add_argument("--phylum", help="phylum for all records (default: from each record's taxonomy)")
    a = ap.parse_args()
    if a.cmd == "build":
        ex = {x.strip().split(".")[0] for x in a.exclude.read_text().split()} if a.exclude else set()
        for r in build(a.gb, a.out, ex, a.cpu):
            print("\t".join(str(r[k]) for k in ("phylum", "gene", "n_ref", "n_used", "median_len")))
        return 0
    rows = search_dir(a.input, a.lib, a.phylum)
    with open(a.input / "rrna_hmm_report.tsv", "w") as h:
        w = csv.DictWriter(h, fieldnames=["file", "record", "phylum", "gene", "location", "status"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{sum(r['status'] == 'added' for r in rows)} rRNA(s) added; report {a.input / 'rrna_hmm_report.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
