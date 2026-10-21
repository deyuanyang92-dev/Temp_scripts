#!/usr/bin/env python3
"""mito_rescue.py - add genes that MitoZ missed, taken from a MITOS2 annotation of the same sequence.

MitoZ gives the most exact protein-gene boundaries, but misses some short or divergent genes that
MITOS2 finds. 2026-10-20 dev benchmark (NCBI NC_ annotation as reference, offline simulation):
  MitoZ-supported phyla (45 records): protein genes found 579 -> 580 (atp8), rRNA 81 -> 87, no new FP.
  Other phyla (48 records, MitoZ with a borrowed protein database): protein genes 525 -> 548, rRNA 75 -> 82.
MITOS2 boundaries are less exact (52 % vs 80 % both ends exact), so nothing MitoZ annotated is changed:
only genes absent from the record are added (each listed in rescue_report.tsv with a hint to check it).

Rules (each from the benchmark and from the 2026-10-21 review):
  - protein gene: added when no CDS of that gene exists and the MITOS2 gene does not lie mostly
    (> 50 %) on an existing CDS, rRNA or tRNA, or on a gene added earlier in the same call (longest first;
    that is a naming disagreement, not a missed gene).
    Containment release: an existing CDS/rRNA stops blocking only when it looks like an intron-spanning or
    merged call: it is longer than 1.5 x MAX_LEN of its own gene (MAX_LEN = longest expected length in nt)
    AND wholly contains >= 2 MITOS2 genes of other names that are absent from the record (NC_006902.1:
    MitoZ nad5 over its 11 kb group I intron, 10 genes inside). A normal-length gene never releases, so
    homologue fragments (nad4 / nad2 inside nad5) are not added as bogus CDS.
  - rRNA: MITOS2 often splits a divergent rRNA into short pieces (29-200 nt); the pieces (-a/-b) are not
    joined, only the longest piece of each rRNA, and only when >= 300 nt, is added (rescued-rRNA FP 13 -> 3
    on the other phyla). CDS fragments (-a/-b) are joined only on the same strand with a gap of 0-3000 nt.
  - /translation of a rescued CDS: trimmed to whole codons, first codon M if it is a start codon of the
    table, trailing stop removed; omitted when an internal stop remains.
  - tRNAs are not handled here (MITOS2 already annotates them in the default pipeline).
  - The feature /note is neutral (NOTE); the advice to check the boundaries is the "hint" column of
    rescue_report.tsv.

MITOS2: Donath et al. 2019, Nucleic Acids Res 47:10543 (doi:10.1093/nar/gkz833). It runs through
scripts/runmitos_fast.py (Infernal searches in parallel, same results as runmitos).

Usage: python3 mito_rescue.py -i 04.final_gb/ --genetic_code 5 [--mitos_bin DIR] [--mitos_refdir DIR]
Writes the GenBank files in place and rescue_report.tsv next to them.
"""
import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402

MITOS_TRNA = re.compile(r"^trn([A-Z])(\d?)(?:_(\d+))?(?:-([a-z]))?\(([a-z-]{3})\)", re.I)
MITOS_GENE = re.compile(r"^(.+?)(?:_(\d+))?(?:-([a-z]))?$")       # nad5, nad5_1 (2nd copy), cox2-a (fragment)
MITOS_SKIP = re.compile(r"^(OH|OL)", re.I)
MIN_RRNA = 300

# MITOS2 name -> MitoZ gene name and product (as in MitoZ's final GenBank)
MITOZ_NAMES = {
    "atp6": ("ATP6", "ATP synthase F0 subunit 6"), "atp8": ("ATP8", "ATP synthase F0 subunit 8"),
    "cox1": ("COX1", "cytochrome c oxidase subunit I"), "cox2": ("COX2", "cytochrome c oxidase subunit II"),
    "cox3": ("COX3", "cytochrome c oxidase subunit III"), "cob": ("CYTB", "cytochrome b"),
    "nad1": ("ND1", "NADH dehydrogenase subunit 1"), "nad2": ("ND2", "NADH dehydrogenase subunit 2"),
    "nad3": ("ND3", "NADH dehydrogenase subunit 3"), "nad4": ("ND4", "NADH dehydrogenase subunit 4"),
    "nad4l": ("ND4L", "NADH dehydrogenase subunit 4L"), "nad5": ("ND5", "NADH dehydrogenase subunit 5"),
    "nad6": ("ND6", "NADH dehydrogenase subunit 6"),
    "rrnS": ("s-rRNA", "12S ribosomal RNA"), "rrnL": ("l-rRNA", "16S ribosomal RNA"),
}
# longest expected gene length (nt), by MITOS2 name; a CDS/rRNA is a possible intron-spanning / merged call
# only when longer than 1.5 x this value
MAX_LEN = {
    "cox1": 1700, "cox2": 800, "cox3": 850, "cob": 1250, "nad1": 1050, "nad2": 1150, "nad3": 420,
    "nad4": 1450, "nad4l": 320, "nad5": 2000, "nad6": 600, "atp6": 750, "atp8": 250, "rrnS": 1100, "rrnL": 1900,
}
CONTAINER_FACTOR = 1.5
MAX_FRAG_GAP = 3000      # CDS fragments (-a/-b) further apart than this are not joined
NOTE = "annotated by similarity to a MITOS2 prediction (not found by the primary annotation)"
HINT = "MITOS2 boundaries are less exact than MitoZ's (52 % vs 80 % both ends exact): check the ends manually"
_MAX_LEN_BY_NAME = {n.lower(): MAX_LEN[g] for g, (n, _) in MITOZ_NAMES.items()}
_MAX_LEN_BY_NAME.update({g.lower(): v for g, v in MAX_LEN.items()})


def max_len(gene: str) -> Optional[int]:
    """Longest expected length (nt) of a gene given by its MitoZ or MITOS2 name; None when unknown."""
    return _MAX_LEN_BY_NAME.get(gene.lower())


def _gap(a: FeatureLocation, b: FeatureLocation) -> int:
    """Nucleotides between two intervals (< 0 when they overlap)."""
    return max(int(b.start) - int(a.end), int(a.start) - int(b.end))


def _chains(frags: List[tuple]) -> List[List[FeatureLocation]]:
    """Fragments (-a, -b ...) in letter order -> chains of parts to join: only parts on the same strand
    with a gap of 0..MAX_FRAG_GAP nt to the previous part are joined; others start a new feature."""
    chains: List[List[FeatureLocation]] = []
    for _, parts in sorted(frags, key=lambda x: x[0]):
        last = chains[-1][-1] if chains else None
        if last is not None and last.strand == parts[0].strand and 0 <= _gap(last, parts[0]) <= MAX_FRAG_GAP:
            chains[-1].extend(parts)
        else:
            chains.append(list(parts))
    return chains


def mitos_to_record(result_dir: Path, seq: Seq, acc: str) -> SeqRecord:
    """GenBank-like record from MITOS2 result.bed (0-based). Copies (_0, _1) stay separate features of the
    same gene; CDS fragments (-a, -b) are joined in letter order when on the same strand and at most
    MAX_FRAG_GAP nt apart (otherwise separate features); rRNA fragments are never joined (rescue then takes
    the longest piece); start > end wraps the origin; OH/OL skipped."""
    n = len(seq)
    groups: Dict[tuple, List[tuple]] = {}
    for i, line in enumerate((result_dir / "result.bed").read_text().splitlines()):
        c = line.split("\t")
        if len(c) < 6 or MITOS_SKIP.match(c[3]):
            continue
        strand = -1 if c[5].strip() == "-" else 1
        m = MITOS_TRNA.match(c[3])
        if m:
            key = ("tRNA", f"trn{m.group(1).upper()}{m.group(2)}", m.group(5).lower().strip("-"), m.group(3) or "")
            frag = m.group(4) or ""
        else:
            g = MITOS_GENE.match(c[3])
            base = g.group(1) if g else c[3]
            kind = "rRNA" if base.lower().startswith("rrn") else "CDS"
            # rRNA pieces are separate features: the line number makes each key unique
            key = (kind, base, "", f"line{i}" if kind == "rRNA" else (g.group(2) or "") if g else "")
            frag = (g.group(3) or "") if g else ""
        s0, e0 = int(c[1]), int(c[2])
        if s0 <= e0:
            parts = [FeatureLocation(s0, e0, strand)]
        else:                                    # circular feature across the origin
            a, b = FeatureLocation(s0, n, strand), FeatureLocation(0, e0, strand)
            parts = [b, a] if strand == -1 else [a, b]
        groups.setdefault(key, []).append((frag, parts))
    r = SeqRecord(seq, id=acc)
    r.annotations["molecule_type"] = "DNA"
    for (kind, name, ac, _), frags in groups.items():
        for parts in _chains(frags):
            loc = parts[0] if len(parts) == 1 else CompoundLocation(parts)
            q = {"gene": [name]}
            if kind == "tRNA":
                q["anticodon"] = [ac]
            r.features.append(SeqFeature(loc, type="gene", qualifiers={"gene": [name]}))
            r.features.append(SeqFeature(loc, type=kind, qualifiers=q))
    return r


def _positions(f: SeqFeature, n: int) -> set:
    return {i % n for p in f.location.parts for i in range(int(p.start), int(p.end))}


def _translation(f: SeqFeature, seq: Any, code: int) -> Optional[str]:
    """Protein of a rescued CDS: whole codons only, first codon M when it is a start codon of the table,
    one trailing stop removed. None when an internal stop remains (the call is doubtful: no /translation)."""
    nt = f.extract(seq)
    nt = nt[: len(nt) - len(nt) % 3]
    if not len(nt):
        return None
    aa = str(nt.translate(table=code))
    if aa.endswith("*"):
        aa = aa[:-1]
    if not aa or "*" in aa:
        return None
    try:
        starts = CodonTable.unambiguous_dna_by_id[code].start_codons
    except KeyError:
        starts = []
    return "M" + aa[1:] if str(nt[:3]).upper() in starts else aa


def rescue(rec: SeqRecord, mitos: SeqRecord, code: int) -> List[str]:
    """Add to rec (MitoZ annotation) the protein genes and rRNAs of mitos (mitos_to_record) that rec lacks.
    Returns the MitoZ-style names added."""
    n = len(rec)
    have = {f.qualifiers.get("gene", [""])[0] for f in rec.features if f.type in ("CDS", "rRNA")}
    others = [(MITOZ_NAMES[f.qualifiers["gene"][0]][0], _positions(f, n)) for f in mitos.features
              if f.type in ("CDS", "rRNA") and f.qualifiers["gene"][0] in MITOZ_NAMES]
    taken = []
    for f in rec.features:
        if f.type not in ("CDS", "rRNA", "tRNA"):
            continue
        pos = _positions(f, n)
        name = f.qualifiers.get("gene", [""])[0]
        limit = max_len(name)
        if f.type != "tRNA" and limit is not None and len(pos) > CONTAINER_FACTOR * limit:
            # far longer than the gene can be: may span an intron or merge genes; releases the genes inside
            # it that are absent from the record, when there are >= 2 of them
            inside = {g for g, p in others if g != name and g not in have and p <= pos}
            if len(inside) >= 2:
                continue
        taken.append(pos)
    best: Dict[str, SeqFeature] = {}
    for f in mitos.features:
        if f.type not in ("CDS", "rRNA") or f.qualifiers["gene"][0] not in MITOZ_NAMES:
            continue
        name = MITOZ_NAMES[f.qualifiers["gene"][0]][0]
        pos = _positions(f, n)
        if name in have or any(len(pos & t) > 0.5 * len(pos) for t in taken):
            continue
        if f.type == "rRNA" and len(pos) < MIN_RRNA:
            continue
        if name not in best or len(pos) > len(_positions(best[name], n)):
            best[name] = f
    added = []
    # longest first; every accepted gene blocks later candidates lying mostly on it
    for name, f in sorted(best.items(), key=lambda kv: -len(_positions(kv[1], n))):
        pos = _positions(f, n)
        if any(len(pos & t) > 0.5 * len(pos) for t in taken):
            continue
        taken.append(pos)
        mname = f.qualifiers["gene"][0]
        product = MITOZ_NAMES[mname][1]
        q = {"gene": [name], "product": [product], "note": [NOTE]}
        if f.type == "CDS":
            q = {"gene": [name], "codon_start": ["1"], "transl_table": [str(code)], "product": [product],
                 "note": [NOTE]}
            tr = _translation(f, rec.seq, code)
            if tr:
                q["translation"] = [tr]
        rec.features += [SeqFeature(f.location, type="gene", qualifiers={"gene": [name]}),
                         SeqFeature(f.location, type=f.type, qualifiers=q)]
        added.append(name)
    if added:
        src = [f for f in rec.features if f.type == "source"]
        rec.features = src + sorted((f for f in rec.features if f.type != "source"),
                                    key=lambda f: (int(f.location.start), f.type != "gene"))
    return added


def run_mitos(fa: Path, code: int, out: Path, linear: bool, mitos_bin: Optional[str] = None,
              refdir: Optional[str] = None, ref: str = "refseq89m") -> Path:
    """MITOS2 protein-gene and rRNA annotation (no tRNA, intron, OH/OL) via scripts/runmitos_fast.py."""
    b = mito_env.mitos_bin(mitos_bin)
    d = mito_env.mitos_refdir(refdir)
    if not (b and d):
        raise RuntimeError("MITOS2 not found (install.sh, or --mitos_bin / --mitos_refdir)")
    out.mkdir(parents=True, exist_ok=True)
    cmd = ["env", f"PATH={b}:/usr/bin:/bin", str(Path(b) / "python"),
           str(Path(__file__).resolve().parent / "scripts" / "runmitos_fast.py"),
           "-i", str(fa), "-c", str(code), "-o", str(out), "-r", ref, "--refdir", d,
           "--trna", "0", "--intron", "0", "--oril", "0", "--orih", "0", "--noplots", *(["--linear"] if linear else [])]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if not (out / "result.bed").exists():
        raise RuntimeError(f"MITOS2 produced no result.bed: {(p.stderr or p.stdout)[-300:]}")
    return out


def rescue_dir(final_dir: Path, code: int, run: Optional[Callable] = None, **kw) -> List[dict]:
    """Rescue every record of every *.gbf in final_dir (in place); returns report rows."""
    if run is None:
        def run(fa, c, out, linear):
            return run_mitos(fa, c, out, linear, **kw)
    rows = []
    for gbf in sorted(Path(final_dir).rglob("*.gbf")):
        recs = list(SeqIO.parse(str(gbf), "genbank"))
        changed = False
        for rec in recs:
            have = {f.qualifiers.get("gene", [""])[0] for f in rec.features if f.type in ("CDS", "rRNA")}
            if all(v[0] in have for v in MITOZ_NAMES.values()):      # nothing to rescue: MITOS2 not run
                continue
            with tempfile.TemporaryDirectory() as td:
                fa = Path(td) / "in.fa"
                fa.write_text(f">x\n{rec.seq}\n")
                (Path(td) / "mitos").mkdir()
                try:
                    res = run(fa, code, Path(td) / "mitos", rec.annotations.get("topology") == "linear")
                except RuntimeError as e:
                    rows.append({"file": gbf.name, "record": rec.id, "gene": "", "location": "", "status": f"failed: {e}",
                                 "hint": ""})
                    continue
                added = rescue(rec, mitos_to_record(res, rec.seq, rec.id), code)
            for name in added:
                f = next(x for x in rec.features if x.type != "gene" and x.qualifiers.get("gene") == [name])
                rows.append({"file": gbf.name, "record": rec.id, "gene": name, "location": str(f.location), "status": "added",
                             "hint": HINT})
            changed = changed or bool(added)
        if changed:
            SeqIO.write(recs, str(gbf), "genbank")
    return rows


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="directory of final GenBank files (*.gbf), edited in place")
    p.add_argument("--genetic_code", type=int, required=True)
    p.add_argument("--mitos_bin")
    p.add_argument("--mitos_refdir")
    a = p.parse_args(argv)
    if not shutil.which("env"):
        print("[mito_rescue] 'env' not found", file=sys.stderr)
        return 2
    rows = rescue_dir(Path(a.input), a.genetic_code, mitos_bin=a.mitos_bin, refdir=a.mitos_refdir)
    rep = Path(a.input) / "rescue_report.tsv"
    with open(rep, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["file", "record", "gene", "location", "status", "hint"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n = sum(r["status"] == "added" for r in rows)
    print(f"[mito_rescue] {n} gene(s) added from MITOS2; {sum(r['status'] != 'added' for r in rows)} failure(s) -> {rep}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
