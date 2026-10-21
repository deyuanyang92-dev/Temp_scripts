#!/usr/bin/env python3
"""mito_rrna_hmm.py - last-resort 12S/16S search with per-phylum profile HMMs (MitoM).

Most genomes are annotated well by the cores (MitoZ / MITOS2 covariance models); this step runs only for
an rRNA that is still missing after them and the MITOS2 rescue. It searches the genome with an nhmmer
profile HMM built from the NCBI rRNAs of the record's own phylum, as EMMA does with its (vertebrate)
HMMs (knowledge/tools/emma.md).
Evidence (2026-10-21 dev benchmark, 186 NCBI rRNAs of 93 records, leave-genus-out, knowledge/topics/
rrna_annotation_methods.md section 6): at search level MITOS2 CM (fragments merged) found 166, CM then
per-phylum HMM 175 (knowledge/experiments/rrna_ref/rr_detect_summary_dev.txt). That experiment merged
all same-strand hits within 300 nt and applied none of the filters below, so it is not a measurement of
this code. In the full pipeline (HMM step before mito_curate, library without any benchmark genus;
pipeline_ab_dev.txt, run with an earlier version of this step) missed rRNAs went 17 -> 14, false
positives 9 -> 10, exactly matching ones 59 -> 59 (present but partial rRNAs count as present). The
coverage / partial-end / window rules below were added after that run and are covered by unit tests only.

Rules:
  - only rRNAs absent from the record are searched (names as in mito_curate.RRNA_NAMES and free text);
    the phylum is taken from the record's taxonomy (or --phylum); phyla with fewer than MIN_REFS
    reference rRNAs have no model and are skipped;
  - hits of one model are joined only when they are pieces of one gene: same strand, disjoint and
    consecutive model columns in the same order as in the genome, <= MAX_GAP nt apart. Hits over the
    same model columns are copies (mitochondrial genes can be duplicated) and are never joined; all such
    clusters are built first and the one with the highest summed score wins. A second copy is only
    reported ("second copy at X..Y not added"), never added;
  - circular records are searched as seq + first 2000 nt, so a gene across the origin is one hit and
    becomes a join location; clusters that start in the appended part are duplicates and dropped;
  - the span must cover >= MIN_COVER (0.4) of the model columns (LENG of the .hmm) and be >= MIN_RRNA
    (200) nt; <= LEN_HI (1.3) x the 95th percentile length of the references (manifest p95_len) and, when
    complete at both ends, >= LEN_LO (0.7) x their 5th percentile; it must not lie mostly (> 50 %) on an
    existing CDS / rRNA. A span whose model columns start after column 30 / end before LENG - 30 is added
    with < / > on the 5' / 3' end ("added (partial 5')"); report columns coverage and note;
  - other ends are not refined here (mito_curate.py R1); the feature /note is a plain statement of the
    method, the "check boundaries manually" hint is only in the report.

Usage:
  python3 scripts/fetch_rrna_refs.py -o refs --all-metazoa [--set benchmark/sets/mito_bench_2026-10.tsv]
  python3 mito_rrna_hmm.py build --gb refs/*.gb -o models/rrna_hmm [--exclude accessions.txt] [--manifest-only]
  python3 mito_rrna_hmm.py search -i 04.final_gb/ [--lib models/rrna_hmm] [--phylum Ctenophora]
build: one GenBank file per phylum (file name = phylum), made by scripts/fetch_rrna_refs.py (reusable form
of knowledge/experiments/rrna_ref/fetch_refs.py, which also writes refs_<phylum>.txt accession lists);
rRNAs of 300-4000 nt with <= 1 % N; --manifest-only recomputes manifest.tsv (lengths, percentiles) without
re-aligning. search writes the GenBank files in place (all *.gbf / *.gb below the directory; failures
become "failed: ..." rows) and rrna_hmm_report.tsv next to them.

HMMER: Eddy 2011, PLoS Comput Biol 7:e1002195 (doi:10.1371/journal.pcbi.1002195); nhmmer: Wheeler &
Eddy 2013, Bioinformatics 29:2487 (doi:10.1093/bioinformatics/btt403); MAFFT: Katoh & Standley 2013,
Mol Biol Evol 30:772 (doi:10.1093/molbev/mst010).
"""
import argparse
import csv
import os
import random
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Callable, Dict, Iterable, List, NamedTuple, Optional, Set, Tuple, Union

from Bio import SeqIO
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, ExactPosition, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402
from mito_curate import RRNA_NAMES, key as name_key  # noqa: E402  (one alias table for rRNA names: srrna, rrn12, rrn16 ...)

LIB = Path(__file__).resolve().parent / "models" / "rrna_hmm"
GENES = {"rrnS": ("s-rRNA", "12S ribosomal RNA"), "rrnL": ("l-rRNA", "16S ribosomal RNA")}
NAMES = {"rrnS": {"rrns", "s-rrna", "12s", "rns", "12s rrna"}, "rrnL": {"rrnl", "l-rrna", "16s", "rnl", "16s rrna"}}
MIN_REFS, MAX_ALN, MIN_RRNA, MAX_GAP, EVALUE, MAX_REF_LEN = 3, 60, 200, 300, 1e-3, 4000
LEN_HI, LEN_LO, LEN_HI_MEDIAN = 1.3, 0.7, 1.5   # span <= LEN_HI x p95 and (if complete) >= LEN_LO x p5 of the references;
                                                  # manifests from before p5/p95 existed: <= LEN_HI_MEDIAN x median
ORIGIN_EXT = 2000   # circular records are searched as seq + seq[:ORIGIN_EXT], so a gene across the origin is contiguous
MANIFEST = ["phylum", "gene", "n_ref", "n_used", "median_len", "p5_len", "p95_len"]
MIN_COVER, PARTIAL_COLS = 0.4, 30   # model fraction a span must cover; columns missing at an end = partial there
NOTE = "annotated by profile HMM similarity (per-phylum rRNA model, {}, {} references)"      # goes into the feature
CHECK = "ends approximate (span of the profile HMM hits, not curated): check boundaries manually"   # report only
ERRORS = (RuntimeError, subprocess.CalledProcessError, OSError, ValueError)   # become "failed: ..." report rows


class Hit(NamedTuple):
    start: int      # 0-based genome start (envelope)
    end: int        # genome end, exclusive
    strand: int
    hmm_from: int   # model columns, 1-based
    hmm_to: int
    score: float


def gene_of(f: SeqFeature) -> Optional[str]:
    """rrnS / rrnL for an rRNA feature named in MitoZ, NCBI or free style; None otherwise. Exact aliases are
    mito_curate.RRNA_NAMES (shared with R3); NAMES / small / large are the looser free-text match."""
    qs = f.qualifiers.get("gene", []) + f.qualifiers.get("product", [])
    alias = next((RRNA_NAMES[name_key(q)] for q in qs if name_key(q) in RRNA_NAMES), None)
    if alias:
        return alias
    s = " ".join(qs).lower()
    for g, keys in NAMES.items():                    # a name must not continue a longer one: "16s-rrna" is not "s-rrna"
        if re.search(r"(?<![a-z0-9])(?:" + "|".join(map(re.escape, keys)) + ")", s) or ("small" in s if g == "rrnS" else "large" in s):
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


def clusters(hits: List[Hit]) -> List[Hit]:
    """Every cluster of hits that are pieces of one gene (see module rules), best summed score first. Each
    unassigned hit, best score first, seeds a cluster that grows by joining the hits that fit it."""
    rest = sorted(hits, key=lambda h: -h.score)
    out = []
    while rest:
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
        out.append(span)
    return sorted(out, key=lambda c: -c.score)


def best_span(hits: List[Hit]) -> Optional[Hit]:
    """The cluster with the highest summed score."""
    return best_spans(hits)[0]


def best_spans(hits: List[Hit], n: Optional[int] = None) -> Tuple[Optional[Hit], Optional[Hit]]:
    """(best cluster, runner-up that is a second copy of it). A copy covers mostly the same model columns
    (>= 50 % overlap) at other genome positions; a fragment over a few columns is not a copy.
    n = genome length of a circular record searched with ORIGIN_EXT nt appended: clusters that start in the
    appended part duplicate the same genes in the first part and are dropped; positions count modulo n."""
    cl = [c for c in clusters(hits) if n is None or c.start < n]
    if not cl:
        return None, None
    best = cl[0]
    used = {i % n if n else i for i in range(best.start, best.end)}

    def copy_of(c: Hit) -> bool:
        both = min(best.hmm_to, c.hmm_to) - max(best.hmm_from, c.hmm_from) + 1
        union = max(best.hmm_to, c.hmm_to) - min(best.hmm_from, c.hmm_from) + 1
        return not used & {i % n if n else i for i in range(c.start, c.end)} and both >= 0.5 * union
    return best, next((c for c in cl[1:] if copy_of(c)), None)


def read_manifest(lib: Path) -> Dict[tuple, dict]:
    f = lib / "manifest.tsv"
    if not f.exists():
        return {}
    with open(f, newline="") as h:
        return {(r["phylum"], r["gene"]): r for r in csv.DictReader(h, delimiter="\t")}


_HEADERS: Dict[tuple, Dict[str, str]] = {}


def _header(hmm: Path) -> Dict[str, str]:
    """Header fields (NAME, LENG, NSEQ ...) of a .hmm file; cached per file version."""
    k = (str(hmm), hmm.stat().st_mtime_ns)
    if k not in _HEADERS:
        h: Dict[str, str] = {}
        with open(hmm) as f:
            for line in f:
                if line.startswith("HMM "):
                    break
                c = line.split(None, 1)
                if len(c) == 2:
                    h[c[0]] = c[1].strip()
        _HEADERS[k] = h
    return _HEADERS[k]


def model_len(hmm: Path) -> int:
    """LENG (number of match columns) from the .hmm header."""
    try:
        return int(_header(hmm)["LENG"])
    except (KeyError, ValueError):
        raise ValueError(f"no LENG line in the header of {hmm}") from None


def _why(e: Exception) -> str:
    """One-line reason for a report row: nhmmer's stderr when the exception carries it."""
    err = getattr(e, "stderr", None)
    if isinstance(err, bytes):
        err = err.decode(errors="replace")
    return " ".join(str(err).split()[-60:] if err else str(e).split())


def nhmmer(hmm: Path, fa: Path) -> List[Hit]:
    exe = mito_env.find("nhmmer")
    if not exe:
        raise RuntimeError("nhmmer not found (HMMER 3; conda env mitotools)")
    with tempfile.TemporaryDirectory() as td:
        tbl = Path(td) / "hits.tbl"
        p = subprocess.run([exe, "--dna", "--cpu", "2", "-E", str(EVALUE), "--tblout", str(tbl), "-o", "/dev/null",
                            str(hmm), str(fa)], capture_output=True, text=True)
        if p.returncode:
            raise RuntimeError(f"nhmmer failed (exit {p.returncode}): {' '.join(p.stderr.split()[-60:])}")
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


def _loc(h: Hit, n: int) -> str:
    """1-based report location; a span across the origin of a circular record reads S..n,1..E."""
    rng = f"{h.start + 1}..{h.end}" if h.end <= n else f"{h.start + 1}..{n},1..{h.end - n}"
    return rng + ("+" if h.strand == 1 else "-")


def _partial(span: Hit, leng: int) -> List[str]:
    """Ends of the gene (5', 3') whose model columns are missing from the span."""
    return [e for e, miss in (("5'", span.hmm_from > PARTIAL_COLS), ("3'", span.hmm_to < leng - PARTIAL_COLS)) if miss]


def _location(span: Hit, partial: List[str], n: int) -> Union[FeatureLocation, CompoundLocation]:
    """Location with < / > on the partial ends (5' = genome start on the plus strand, genome end on minus).
    A span across the origin (end > n) becomes a join in transcript order, as mito_rescue writes it."""
    left, right = ("5'", "3'") if span.strand == 1 else ("3'", "5'")          # which gene end lies at the genome start / end
    s = BeforePosition(span.start) if left in partial else ExactPosition(span.start)
    e = AfterPosition(span.end) if right in partial else ExactPosition(span.end)
    if span.end <= n:
        return FeatureLocation(s, e, span.strand)
    e = AfterPosition(span.end - n) if right in partial else ExactPosition(span.end - n)
    first, last = FeatureLocation(s, ExactPosition(n), span.strand), FeatureLocation(ExactPosition(0), e, span.strand)
    return CompoundLocation([first, last] if span.strand == 1 else [last, first])


def _verdict(span: Hit, m: dict, leng: int, taken: List[Set[int]], n: int) -> str:
    """"added" / "added (partial 5')" ... or the reason the joined span is rejected."""
    length, pos = span.end - span.start, {i % n for i in range(span.start, span.end)}
    cover, part = (span.hmm_to - span.hmm_from + 1) / leng, _partial(span, leng)
    hi = LEN_HI * float(m["p95_len"]) if m.get("p95_len") else LEN_HI_MEDIAN * float(m["median_len"])
    lo = LEN_LO * float(m["p5_len"]) if m.get("p5_len") and not part else 0   # a partial gene is short by definition
    if cover < MIN_COVER:
        return f"rejected: low coverage ({cover:.2f} of the model < {MIN_COVER})"
    if length < max(MIN_RRNA, lo):
        return f"rejected: short ({length} nt < {max(MIN_RRNA, lo):.0f})"
    if length > hi:
        return f"rejected: long ({length} nt > {hi:.0f} = {LEN_HI if m.get('p95_len') else LEN_HI_MEDIAN} x reference p95/median)"
    if any(len(pos & t) > 0.5 * length for t in taken):
        return "rejected: lies on an existing gene"
    return f"added (partial {' and '.join(part)})" if part else "added"


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
    rows: List[dict] = []
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "q.fa"
        circ = rec.annotations.get("topology") == "circular"
        seq = str(rec.seq)
        fa.write_text(f">{rec.id}\n{seq + seq[:min(ORIGIN_EXT, n)] if circ else seq}\n")
        for g in want:
            name, product = GENES[g]
            m = man.get((phy, g))
            if not m:
                rows.append(dict(row, gene=name, location="", status="skipped: no model for this gene"))
                continue
            hmm = lib / f"{phy}.{g}.hmm"
            try:
                span, copy = best_spans(run(hmm, fa), n if circ else None)
                leng = model_len(hmm)
            except ERRORS as e:
                rows.append(dict(row, gene=name, location="", status=f"failed: {_why(e)}"))
                continue
            if span is None:
                rows.append(dict(row, gene=name, location="", status="not found"))
                continue
            status = _verdict(span, m, leng, taken, n)
            rows.append(dict(row, gene=name, location=_loc(span, n), status=status, note=CHECK if status.startswith("added") else "",
                             coverage=f"{(span.hmm_to - span.hmm_from + 1) / leng:.2f}"))
            if not status.startswith("added"):
                continue
            fl = _location(span, _partial(span, leng), n)
            q = {"gene": [name], "product": [product], "note": [NOTE.format(phy, m["n_used"])]}
            rec.features += [SeqFeature(fl, type="gene", qualifiers={"gene": [name]}), SeqFeature(fl, type="rRNA", qualifiers=q)]
            taken.append({i % n for i in range(span.start, span.end)})
            if copy and _verdict(copy, m, leng, taken, n).startswith("added"):
                rows.append(dict(row, gene=name, location=_loc(copy, n), status=f"second copy at {_loc(copy, n)[:-1]} not added"))
    if any(r["status"].startswith("added") for r in rows):
        src = [f for f in rec.features if f.type == "source"]
        rec.features = src + sorted((f for f in rec.features if f.type != "source"),
                                    key=lambda f: (int(f.location.start), f.type != "gene"))
    return rows


def _write_gb(recs: List[SeqRecord], path: Path) -> None:
    """Write a GenBank file atomically: temp file in the same directory, then os.replace()."""
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as h:
            SeqIO.write(recs, h, "genbank")
        shutil.copymode(path, tmp)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def search_dir(final_dir: Path, lib: Path = LIB, phylum: Optional[str] = None,
               run: Optional[Callable[[Path, Path], List[Hit]]] = None) -> List[dict]:
    """All *.gbf / *.gb under final_dir (subdirectories too, e.g. no_cox1_gene/), rewritten in place when an
    rRNA was added. Never raises for one bad record or file: that is a "failed: ..." row."""
    rows = []
    for f in sorted(set(final_dir.rglob("*.gbf")) | set(final_dir.rglob("*.gb"))):
        name, none = f.relative_to(final_dir).as_posix(), {"record": "", "phylum": "", "gene": "", "location": ""}
        try:
            recs = list(SeqIO.parse(str(f), "genbank"))
        except ERRORS as e:
            rows.append(dict(none, file=name, status=f"failed: cannot read ({_why(e)})"))
            continue
        frows: List[dict] = []
        for r in recs:
            try:
                rs = search_record(r, lib, phylum, run)
            except ERRORS as e:
                rs = [dict(none, record=r.id, status=f"failed: {_why(e)}")]
            frows += [dict(x, file=name) for x in rs]
        if any(x["status"].startswith("added") for x in frows):
            try:
                _write_gb(recs, f)
            except ERRORS as e:
                for x in frows:
                    if x["status"].startswith("added"):
                        x["status"] = f"failed: not saved ({_why(e)})"
        rows += frows
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
            if 300 <= len(s) <= MAX_REF_LEN and s.count("N") <= 0.01 * len(s):
                out[g].append(s)
    return out


def _pct(xs: List[int], q: float) -> int:
    """q-quantile (0..1) by linear interpolation between order statistics (numpy default)."""
    v = sorted(xs)
    k = (len(v) - 1) * q
    i = int(k)
    return round(v[i] + (v[min(i + 1, len(v) - 1)] - v[i]) * (k - i))


def build(gbs: List[Path], out: Path, exclude: Optional[Set[str]] = None, cpu: int = 4,
          manifest_only: bool = False) -> List[dict]:
    """One HMM per phylum (GenBank file stem) and rRNA: MAFFT alignment of <= MAX_ALN random references,
    hmmbuild --dna. Writes out/<phylum>.<rrnS|rrnL>.hmm and out/manifest.tsv. manifest_only: recompute
    manifest.tsv from the references only, keeping the existing .hmm files (no mafft / hmmbuild)."""
    mafft, hmmbuild = mito_env.find("mafft"), mito_env.find("hmmbuild")
    if not (manifest_only or (mafft and hmmbuild)):
        raise RuntimeError("mafft and hmmbuild are needed to build the library")
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for gb in sorted(gbs):
        for g, seqs in _ref_rrnas(gb, exclude or set()).items():
            if len(seqs) < MIN_REFS:
                continue
            hmm = out / f"{gb.stem}.{g}.hmm"
            if manifest_only:
                if not hmm.exists():
                    continue
                n_used = int(_header(hmm).get("NSEQ", min(len(seqs), MAX_ALN)))
            else:
                sub = seqs if len(seqs) <= MAX_ALN else random.Random(1).sample(seqs, MAX_ALN)
                n_used = len(sub)
                with tempfile.TemporaryDirectory() as td:
                    fa, aln = Path(td) / "in.fa", Path(td) / "aln.fa"
                    fa.write_text("".join(f">r{i}\n{s}\n" for i, s in enumerate(sub)))
                    aln.write_text(subprocess.run([str(mafft), "--auto", "--thread", str(cpu), "--quiet", str(fa)],
                                                  capture_output=True, text=True, check=True).stdout)
                    subprocess.run([str(hmmbuild), "--dna", "--cpu", str(cpu), "-n", f"{gb.stem}.{g}", str(hmm), str(aln)],
                                   capture_output=True, check=True)
            lens = [len(s) for s in seqs]
            rows.append({"phylum": gb.stem, "gene": g, "n_ref": len(seqs), "n_used": n_used,
                         "median_len": int(statistics.median(lens)), "p5_len": _pct(lens, 0.05), "p95_len": _pct(lens, 0.95)})
    old = read_manifest(out)
    old.update({(r["phylum"], r["gene"]): r for r in rows})
    tmp = out / "manifest.tsv.tmp"
    with open(tmp, "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=MANIFEST, delimiter="\t")
        w.writeheader()
        w.writerows(sorted(old.values(), key=lambda r: (r["phylum"], r["gene"])))
    os.replace(tmp, out / "manifest.tsv")
    return rows


REPORT = ["file", "record", "phylum", "gene", "location", "status", "coverage", "note"]


def write_report(rows: List[dict], path: Path) -> None:
    with open(path, "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=REPORT, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="build per-phylum rRNA HMMs from reference GenBank files")
    b.add_argument("--gb", nargs="+", type=Path, required=True, help="one GenBank file per phylum (stem = phylum)")
    b.add_argument("-o", "--out", type=Path, default=LIB)
    b.add_argument("--exclude", type=Path, help="file with accessions to leave out (one per line)")
    b.add_argument("--cpu", type=int, default=4)
    b.add_argument("--manifest-only", action="store_true",
                   help="recompute manifest.tsv from the references; keep the existing .hmm files (no alignment)")
    s = sub.add_parser("search", help="search missing rRNAs in final GenBank files (in place)")
    s.add_argument("-i", "--input", type=Path, required=True, help="directory of *.gbf / *.gb")
    s.add_argument("--lib", type=Path, default=LIB)
    s.add_argument("--phylum", help="phylum for all records (default: from each record's taxonomy)")
    a = ap.parse_args()
    if a.cmd == "build":
        ex = {x.strip().split(".")[0] for x in a.exclude.read_text().split()} if a.exclude else set()
        for r in build(a.gb, a.out, ex, a.cpu, a.manifest_only):
            print("\t".join(str(r[k]) for k in MANIFEST))
        return 0
    rows = search_dir(a.input, a.lib, a.phylum)
    write_report(rows, a.input / "rrna_hmm_report.tsv")
    print(f"{sum(r['status'].startswith('added') for r in rows)} rRNA(s) added; report {a.input / 'rrna_hmm_report.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
