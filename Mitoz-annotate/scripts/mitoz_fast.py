#!/usr/bin/env python3
"""Run MitoZ with its tRNA step done by MITOS2 (fast mode) and without the circos plot.

MitoZ annotate spends ~80 % of its time in MiTFi (mitoz/annotate/annotation.py tRNA_annotation():
`java -jar mitfi.jar -cores 1`, single core because "MiTFi will get no results when use multiple
cores"): 137 of 168 s for Platynereis dumerilii NC_000931. MITOS2 searches the same MiTFi Metazoa
covariance models; in fast mode (scripts/runmitos_fast.py) it takes ~10 s. Here MitoZ's
tRNA_annotation() is replaced by one that runs MITOS2 and writes MiTFi's own result format
(<mt_file>.trna), which MitoZ's own tRNA_ft_mitfi.py turns into the feature table, so control-region
detection, tbl2asn and everything downstream run unchanged. The circos plot (visualize.main, ~8 s)
is skipped. If MITOS2 is not configured or fails, the original MiTFi step runs, so the result is never
missing tRNAs because of this wrapper.

Settings (environment variables, set by batch_mitoz.py):
  YZZ_TRNA          mitos2 (default) | mitoz (original MiTFi, nothing replaced)
  YZZ_MITOS_BIN     MITOS2 environment bin (contains runmitos and python)
  YZZ_MITOS_REFDIR  MITOS2 reference directory (contains refseq89m)
  YZZ_MITOS_REF     reference name (default refseq89m)
  YZZ_CIRCOS        yes = keep the circos plot (default no)
  YZZ_MITOS_MODE    sensitive (default: MITOS2's sensitive pass for missing tRNAs, run in parallel - same
                    tRNAs as runmitos; fast mode missed 28/28 such tRNAs in the 2026-10-20 dev benchmark)
                    | fast (skip that pass) | full (unmodified runmitos)
  YZZ_PHASE1_FAST   1 = fast tRNA mode in batch_mitoz round 1 (01.mitoz_anno1), which only orients the
                    sequence when round 2 re-annotates it
  YZZ_CLADE_DB      protein database MitoZ ships but does not offer as --clade (Cnidaria,
                    Platyhelminthes-flatworms, Chaetognatha-arrow-worms, Animal): MitoZ runs with
                    --clade Arthropoda and a --profiles_dir that points that name at this database
  YZZ_TRNA_WRAP     file collecting MITOS2 tRNAs across the sequence origin (MiTFi's table cannot hold
                    them); batch_mitoz.py adds them to the final GenBank (mito_trna.add_wrapped)

Usage: <MitoZ env>/bin/python mitoz_fast.py annotate <mitoz annotate arguments>
"""
import os
import re
import subprocess
import sys
import tempfile
import types
from pathlib import Path

from Bio import SeqIO
from mitoz.annotate import annotation
from mitoz.MitoZ import main

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import mito_trna  # noqa: E402  (same directory as batch_mitoz.py)

ORIGINAL_TRNA = annotation.tRNA_annotation
TRNA = os.environ.get("YZZ_TRNA", "mitos2").lower()
MITOS_BIN = os.environ.get("YZZ_MITOS_BIN", "")
MITOS_REFDIR = os.environ.get("YZZ_MITOS_REFDIR", "")
MITOS_REF = os.environ.get("YZZ_MITOS_REF", "refseq89m")
WRAP_FILE = os.environ.get("YZZ_TRNA_WRAP", "")
MITOS_MODE = os.environ.get("YZZ_MITOS_MODE", "sensitive")
CLADE_DB = os.environ.get("YZZ_CLADE_DB", "")


def mitos2_ready() -> bool:
    return bool(MITOS_BIN and MITOS_REFDIR and (Path(MITOS_BIN) / "runmitos").exists()
                and (Path(MITOS_REFDIR) / MITOS_REF).is_dir())


def mitos2_trna_annotation(workdir=None, mt_file=None, genetic_code=None, logger=None):
    """Drop-in for annotation.tRNA_annotation(): same arguments, same output files."""
    workdir = os.path.abspath(workdir)
    os.makedirs(workdir, exist_ok=True)
    mt_file_base = os.path.basename(mt_file)
    trna_out = os.path.join(workdir, mt_file_base + ".trna")
    opts = types.SimpleNamespace(mitos_mode=MITOS_MODE, mitos_ref=MITOS_REF, mitos_refdir=MITOS_REFDIR)
    rows = []
    try:
        for rec in SeqIO.parse(mt_file, "fasta"):
            with tempfile.TemporaryDirectory() as td:
                fa = Path(td) / "g.fa"
                fa.write_text(f">{rec.id}\n{rec.seq}\n")
                preds = mito_trna.mitos2(fa, int(genetic_code), {"runmitos": str(Path(MITOS_BIN) / "runmitos")},
                                         Path(td), opts)
            preds, wrapped = mito_trna.split_wrapped(preds, str(rec.seq))
            for w in wrapped:                     # MiTFi format cannot express a tRNA across the origin
                logger.warning(f"[mitoz_fast] {rec.id}: tRNA {w['aa']} across the origin; "
                               f"added to the final GenBank afterwards ({WRAP_FILE or 'YZZ_TRNA_WRAP not set: lost'})")
                if WRAP_FILE:
                    with open(WRAP_FILE, "a") as fh:
                        print(rec.id, w["aa"], w["anticodon"], w["seq"], sep="\t", file=fh)
            for p in sorted(preds, key=lambda x: x.start):
                ac = p.anticodon.upper().replace("T", "U") or "---"      # MiTFi writes --- for an unreadable anticodon
                # MiTFi columns: seqid start stop score evalue AC AA model strand (start < stop on both strands)
                rows.append([rec.id, p.start, p.end, p.score or "NA", "NA", ac, p.aa,
                             f"Metazoa_{p.aa}.cm", "+" if p.strand == 1 else "-"])
    except RuntimeError as exc:
        logger.warning(f"[mitoz_fast] MITOS2 failed ({exc}); falling back to MiTFi")
        return ORIGINAL_TRNA(workdir=workdir, mt_file=mt_file, genetic_code=genetic_code, logger=logger)
    with open(trna_out, "w") as fh:
        print("#header\tstart\tstop\tscore\tevalue\tAC\tAA\tmodel\tstrand", file=fh)
        for r in rows:
            print(*r, sep="\t", file=fh)
    logger.info(f"[mitoz_fast] MITOS2 fast: {len(rows)} tRNA(s) -> {trna_out}")
    trna_ft = trna_out + ".ft"
    subprocess.check_call([sys.executable, os.path.join(annotation.annotate_script_dir, "tRNA_ft_mitfi.py"),
                           trna_out, trna_ft])
    return trna_ft


def no_circos(args):
    os.makedirs(args.outdir, exist_ok=True)
    return args.outdir


def gather_existing(*files, logger=None, result_wdir=None):
    """gather_result() without the 'cp' error for the circos files that were not drawn."""
    keep = [f for f in files if f and (os.path.exists(f) or any(c in f for c in "*?["))]
    if keep:
        ORIGINAL_GATHER(*keep, logger=logger, result_wdir=result_wdir)


ORIGINAL_GATHER = annotation.gather_result

def use_clade_db(argv):
    """--clade X (X in YZZ_CLADE_DB) -> --clade Arthropoda --profiles_dir <workdir>/mitom_profiles/X."""
    import mito_env
    work = Path(argv[argv.index("--workdir") + 1]) if "--workdir" in argv else Path.cwd()
    src = Path(annotation.__file__).resolve().parents[1] / "profiles"
    d = mito_env.mitoz_profiles(src, CLADE_DB, work / "mitom_profiles", as_clade="Arthropoda")
    if "--clade" in argv:
        argv[argv.index("--clade") + 1] = "Arthropoda"
    else:
        argv += ["--clade", "Arthropoda"]
    argv += ["--profiles_dir", str(d)]


if __name__ == "__main__":
    # round 1 of batch_mitoz only orients the sequence (its tRNAs are discarded when round 2 runs): fast tRNAs
    wd = sys.argv[sys.argv.index("--workdir") + 1] if "--workdir" in sys.argv else ""
    if os.environ.get("YZZ_PHASE1_FAST") == "1" and "01.mitoz_anno1" in wd:
        MITOS_MODE = "fast"
    if CLADE_DB:
        use_clade_db(sys.argv)
    if TRNA == "mitos2":
        if mitos2_ready():
            annotation.tRNA_annotation = mitos2_trna_annotation
        else:
            print("[mitoz_fast] MITOS2 not configured (YZZ_MITOS_BIN / YZZ_MITOS_REFDIR); using MiTFi",
                  file=sys.stderr)
    if os.environ.get("YZZ_CIRCOS", "no").lower() != "yes":
        annotation.visualize = types.SimpleNamespace(main=no_circos)
        annotation.gather_result = gather_existing
    sys.argv[0] = re.sub(r"(-script\.pyw?|\.exe)?$", "", "mitoz")
    sys.exit(main())
