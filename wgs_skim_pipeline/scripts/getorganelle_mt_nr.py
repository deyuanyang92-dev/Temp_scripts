#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
meangs2get_snakemake_v4.py
==========================
Snakemake scheduler wrapper for batch organelle / rRNA assembly FROM READS.

Author  : (your lab)
Version : 4.0
License : MIT

────────────────────────────────────────────────────────────────────────────
MODES (subcommands)
────────────────────────────────────────────────────────────────────────────
  mgs2get      MEANGS (per-sample) ➜ GetOrganelle (MT) ➜ check/summary
  getorganelle GetOrganelle (MT only)  ➜ check/summary
  nr-only      GetOrganelle (NR/rRNA, -F anonym + -s + --genes) ➜ check/summary
  mt-nr        GetOrganelle (MT) + GetOrganelle (NR/rRNA) ➜ check/summary  [RECOMMENDED]
  all          MEANGS ➜ GetOrganelle (MT) ➜ GetOrganelle (NR/rRNA) ➜ check/summary

────────────────────────────────────────────────────────────────────────────
MULTI-SUFFIX READS DETECTION  (new in v4)
────────────────────────────────────────────────────────────────────────────
  --suffix_fq accepts multiple R1,R2 suffix pairs separated by semicolons
  or as repeated arguments.  Examples:

    # single pair (default already includes common patterns)
    --suffix_fq "_1.clean.fq.gz,_2.clean.fq.gz"

    # multiple pairs (semicolon-separated in one string)
    --suffix_fq "_1.clean.fq.gz,_2.clean.fq.gz;_1.fq.gz,_2.fq.gz;_R1.fastq.gz,_R2.fastq.gz"

    # multiple pairs (space-separated on command line)
    --suffix_fq "_1.clean.fq.gz,_2.clean.fq.gz" "_1.fq.gz,_2.fq.gz"

  Pairs are tried IN ORDER; first match wins per sample name.
  Default pairs:
    1) _1.clean.fq.gz  ,  _2.clean.fq.gz
    2) _1.fq.gz        ,  _2.fq.gz
    3) _R1.fastq.gz    ,  _R2.fastq.gz
    4) _R1.fq.gz       ,  _R2.fq.gz
    5) _1.fastq.gz     ,  _2.fastq.gz

────────────────────────────────────────────────────────────────────────────
GETORGANELLE SPEED TIPS (new in v4)
────────────────────────────────────────────────────────────────────────────
  MT assembly is dramatically faster when you provide a seed sequence (-s):
    --mt_seed_fasta /path/to/related_species_mt.fasta

  Other speed knobs:
    --w_mt 105          word size for MT graph (larger = faster, less sensitive)
    --max_reads_mt 1500000   cap reads used (1.5 M often sufficient for MT)
    --R_mt 5            fewer extension rounds (10 is conservative)

  NR / rRNA assembly REQUIRES both -s and --genes (see GetOrganelle wiki):
    --nr_seed_fasta /path/to/seed.fasta
    --nr_genes_fasta /path/to/label.fasta
    --P_nr 0            RECOMMENDED for short high-copy loci (18S/28S)
    --max_extending_len_nr 250   limit extension length for rRNA

────────────────────────────────────────────────────────────────────────────
ROBUSTNESS DESIGN
────────────────────────────────────────────────────────────────────────────
  - Per-sample jobs NEVER crash the whole workflow:
    Each job writes a stamp file (.mt.done / .nr.done / .meangs.done) with
    ok=0/1 and reason.  Snakemake continues; failures are aggregated into:
      <out_root>/mt_failed.txt
      <out_root>/nr_failed.txt
      <out_root>/failed_all.txt

  - Circular detection: a FASTA is circular if ANY header contains "circular"
    (case-insensitive).

────────────────────────────────────────────────────────────────────────────
SUMMARY OUTPUTS
────────────────────────────────────────────────────────────────────────────
    <out_root>/summary_fasta_mt.txt    normalised MT FASTA (all samples)
    <out_root>/summary_fasta_nr.txt    normalised NR FASTA (all samples)
    <out_root>/summary_fasta_all.txt   MT + NR concatenated
    <out_root>/summary_report.tsv      per-sample/locus stats (TSV)

────────────────────────────────────────────────────────────────────────────
QUICK START EXAMPLES
────────────────────────────────────────────────────────────────────────────

  # 0) Preview (dry-run)
  python3 meangs2get_snakemake_v4.py mt-nr -i 01.cleandata -o out --dryrun

  # 1) MT + NR  (most common, recommended)
  python3 meangs2get_snakemake_v4.py mt-nr \\
      -i 01.cleandata -o out \\
      --nr_seed_fasta /db/rRNA/seed.fasta \\
      --nr_genes_fasta /db/rRNA/label.fasta

  # 2) MT only, with a seed from a related species (faster)
  python3 meangs2get_snakemake_v4.py getorganelle \\
      -i 01.cleandata -o out \\
      --mt_seed_fasta /db/Mytilus_mt.fasta \\
      --max_reads_mt 1500000

  # 3) NR / rRNA only
  python3 meangs2get_snakemake_v4.py nr-only \\
      -i 01.cleandata -o out \\
      --nr_seed_fasta /db/rRNA/seed.fasta \\
      --nr_genes_fasta /db/rRNA/label.fasta

  # 4) Custom read suffixes
  python3 meangs2get_snakemake_v4.py mt-nr \\
      -i 01.cleandata -o out \\
      --suffix_fq "_1.clean.fq.gz,_2.clean.fq.gz;_1.fq.gz,_2.fq.gz" \\
      --nr_seed_fasta /db/seed.fasta \\
      --nr_genes_fasta /db/label.fasta

  # 5) Full MEANGS -> MT -> NR pipeline
  python3 meangs2get_snakemake_v4.py all \\
      -i 01.cleandata -o out \\
      --nr_seed_fasta /db/seed.fasta \\
      --nr_genes_fasta /db/label.fasta \\
      --max_parallel 4
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
MEANGS_EXE_NAME        = "meangs.py"
GETORG_READS_EXE_NAME  = "get_organelle_from_reads.py"

DEFAULT_SUFFIX_PAIRS: List[str] = [
    "_1.clean.fq.gz,_2.clean.fq.gz",
    "_1.fq.gz,_2.fq.gz",
    "_R1.fastq.gz,_R2.fastq.gz",
    "_R1.fq.gz,_R2.fq.gz",
    "_1.fastq.gz,_2.fastq.gz",
]

# GetOrganelle intermediate dirs / files to optionally remove
_GETORG_CLEAN_DIRS  = ["extended_spades", "seed", "tmp"]
_GETORG_CLEAN_FILES = [
    "extended_1_paired.fq",  "extended_2_paired.fq",
    "extended_1_unpaired.fq","extended_2_unpaired.fq",
]


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
def setup_logger(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    for h in list(root.handlers):
        root.removeHandler(h)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(str(log_file), encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def tail_text(s: str, n_lines: int = 120) -> str:
    lines = (s or "").splitlines()
    return "\n".join(lines[-n_lines:])


def safe_name(s: str) -> str:
    s = str(s).strip().replace(" ", "_")
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", s)
    return s or "NA"


def uniq_keep_order(items: Sequence[str]) -> List[str]:
    seen: set = set()
    out: List[str] = []
    for x in items:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def read_list_maybe_txt(raw_vals: Optional[List[str]], base_dir: Path) -> Optional[List[str]]:
    """
    Accept sample-list as:
      --subdir_list A B
      --subdir_list A,B
      --subdir_list list.txt   (one per line; blank / # lines ignored)
    """
    if not raw_vals:
        return None

    if len(raw_vals) == 1:
        p = Path(raw_vals[0])
        if not p.is_file():
            p = base_dir / raw_vals[0]
        if p.is_file():
            items: List[str] = []
            for line in p.read_text(encoding="utf-8", errors="ignore").splitlines():
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                items.append(s)
            return uniq_keep_order(items)

    items = []
    for v in raw_vals:
        for x in str(v).split(","):
            x = x.strip()
            if x:
                items.append(x)
    return uniq_keep_order(items) if items else None


def extra_args_from_list(lst: Optional[List[str]]) -> List[str]:
    out: List[str] = []
    for s in lst or []:
        if s is None:
            continue
        out.extend(shlex.split(str(s)))
    return out


def resolve_path_maybe_relative(p: Optional[str], base: Path) -> Optional[str]:
    if not p:
        return None
    pp = Path(p).expanduser()
    if pp.is_file():
        return str(pp.resolve())
    pp2 = (base / p).expanduser()
    if pp2.is_file():
        return str(pp2.resolve())
    return str(pp.resolve())


# ---------------------------------------------------------------------------
# Multi-suffix parsing  (NEW in v4)
# ---------------------------------------------------------------------------
def parse_suffix_pairs(raw_vals: List[str]) -> List[Tuple[str, str]]:
    """
    Parse paired-end read suffix specifications from command-line arguments.

    Supported formats
    -----------------
    Single pair in one string:
        "_1.clean.fq.gz,_2.clean.fq.gz"

    Multiple pairs via semicolons (one argument):
        "_1.clean.fq.gz,_2.clean.fq.gz;_1.fq.gz,_2.fq.gz"

    Multiple pairs as separate arguments (nargs="+"):
        ["_1.clean.fq.gz,_2.clean.fq.gz", "_1.fq.gz,_2.fq.gz"]

    Mixed:
        ["_1.clean.fq.gz,_2.clean.fq.gz;_1.fq.gz,_2.fq.gz", "_R1.fastq.gz,_R2.fastq.gz"]

    Returns
    -------
    List[Tuple[str, str]]
        Ordered, deduplicated list of (r1_suffix, r2_suffix) tuples.

    Raises
    ------
    ValueError
        If a chunk cannot be split into exactly two non-empty parts.
    """
    pairs: List[Tuple[str, str]] = []
    seen: set = set()

    for item in raw_vals:
        # Split by semicolons or whitespace (but NOT commas, which separate R1/R2)
        chunks = re.split(r"[;\s]+", item.strip())
        for chunk in chunks:
            chunk = chunk.strip()
            if not chunk:
                continue
            parts = chunk.split(",", 1)
            if len(parts) != 2:
                raise ValueError(
                    f"Invalid suffix pair specification: '{chunk}'\n"
                    f"Expected format: R1_SUFFIX,R2_SUFFIX  "
                    f"(e.g. '_1.fq.gz,_2.fq.gz')"
                )
            r1, r2 = parts[0].strip(), parts[1].strip()
            if not r1 or not r2:
                raise ValueError(
                    f"Empty R1 or R2 suffix in: '{chunk}'"
                )
            key = (r1, r2)
            if key not in seen:
                seen.add(key)
                pairs.append(key)

    if not pairs:
        raise ValueError("No valid suffix pairs parsed from input.")
    return pairs


# ---------------------------------------------------------------------------
# Tool resolution
# ---------------------------------------------------------------------------
def which_or_none(exe: str) -> Optional[str]:
    try:
        return shutil.which(exe)
    except Exception:
        return None


def resolve_exe_from_path(exe_name: str, tool_path: Optional[str]) -> Optional[Path]:
    """
    tool_path may be:
      - None
      - a directory  (will look for exe_name inside it)
      - a file path  (the executable itself)
    """
    if not tool_path:
        return None
    p = Path(tool_path).expanduser()
    if p.is_dir():
        return (p / exe_name).resolve()
    return p.resolve()


def build_tool_cmd(
    *,
    tool_label: str,
    exe_name: str,
    base_args: List[str],
    conda_env: Optional[str],
    tool_path: Optional[str],
    conda_exe: str = "conda",
) -> List[str]:
    """
    Build a command list for a tool, respecting:
      1. conda_env   → ``conda run -n ENV exe_name ...``
      2. tool_path   → ``<tool_path>/exe_name ...`` (or ``sys.executable exe_name ...`` for .py)
      3. PATH        → ``shutil.which(exe_name) ...``
    """
    if conda_env:
        conda_bin = which_or_none(conda_exe)
        if not conda_bin:
            raise FileNotFoundError(
                f"[{tool_label}] '{conda_exe}' not found in PATH.\n"
                f"Install conda/mamba or provide the full path via --conda_exe."
            )
        return [conda_bin, "run", "-n", conda_env, "--no-capture-output", exe_name] + base_args

    exe_path = resolve_exe_from_path(exe_name, tool_path)
    if exe_path is not None:
        if not exe_path.exists():
            raise FileNotFoundError(
                f"[{tool_label}] executable not found: {exe_path}\n"
                f"Check --{tool_label.lower()}_path."
            )
        if exe_path.suffix.lower() == ".py":
            return [sys.executable, str(exe_path)] + base_args
        return [str(exe_path)] + base_args

    found = which_or_none(exe_name)
    if not found:
        raise FileNotFoundError(
            f"[{tool_label}] cannot find '{exe_name}' in PATH.\n"
            f"Fixes:\n"
            f"  • conda activate <env> before running this script, OR\n"
            f"  • --conda_env_{tool_label.lower()} <ENV_NAME>, OR\n"
            f"  • --{tool_label.lower()}_path <DIR_OR_SCRIPT>"
        )
    return [found] + base_args


# ---------------------------------------------------------------------------
# Input scanning  (updated for multi-suffix)
# ---------------------------------------------------------------------------
def detect_pairs_in_dir(
    dir_path: Path,
    suffix_pairs: List[Tuple[str, str]],
) -> Dict[str, Tuple[Path, Path]]:
    """
    Scan *dir_path* for paired FASTQ files.

    Each suffix pair (r1_suffix, r2_suffix) is tried in order.
    The first matching suffix pair that produces a valid R1+R2 pair for a
    given sample name wins; later pairs are skipped for that sample name.

    Parameters
    ----------
    dir_path     : directory to scan
    suffix_pairs : ordered list of (r1_suffix, r2_suffix) tuples

    Returns
    -------
    Dict mapping sample_name -> (r1_path, r2_path)
    """
    result: Dict[str, Tuple[Path, Path]] = {}
    if not dir_path.exists():
        return result

    all_files = {f.name: f for f in dir_path.iterdir() if f.is_file()}

    for r1_suffix, r2_suffix in suffix_pairs:
        for fname, fpath in sorted(all_files.items()):
            if not fname.endswith(r1_suffix):
                continue
            stem = fname[: -len(r1_suffix)]
            if not stem:
                continue
            if stem in result:          # already matched by a higher-priority pair
                continue
            r2_name = f"{stem}{r2_suffix}"
            if r2_name in all_files:
                result[stem] = (fpath.resolve(), all_files[r2_name].resolve())

    return result


def collect_tasks(
    input_dir: Path,
    fq_position: str,
    subdir_list: Optional[List[str]],
    suffix_pairs: List[Tuple[str, str]],
) -> List[Tuple[str, Optional[Path], Optional[Path]]]:
    """
    Build the full task list: list of (task_id, r1_path, r2_path).
    r1/r2 are None when reads are missing.
    """
    tasks: List[Tuple[str, Optional[Path], Optional[Path]]] = []

    if fq_position == "flat":
        pairs = detect_pairs_in_dir(input_dir, suffix_pairs)
        allow = set(subdir_list) if subdir_list else None
        for sample, (fq1, fq2) in sorted(pairs.items()):
            if allow and sample not in allow:
                continue
            tasks.append((sample, fq1, fq2))
        return tasks

    # subdir layout
    subdirs = (
        subdir_list
        if subdir_list
        else sorted([p.name for p in input_dir.iterdir() if p.is_dir()])
    )
    for g in subdirs:
        gdir = input_dir / g
        if not gdir.exists() or not gdir.is_dir():
            logging.warning(f"[SKIP] subdir not found: {gdir}")
            continue
        pairs = detect_pairs_in_dir(gdir, suffix_pairs)
        if not pairs:
            tasks.append((g, None, None))
        elif len(pairs) == 1:
            (_sample, (fq1, fq2)) = next(iter(pairs.items()))
            tasks.append((g, fq1, fq2))
        else:
            for sample, (fq1, fq2) in sorted(pairs.items()):
                tasks.append((f"{g}__{sample}", fq1, fq2))
    return tasks


# ---------------------------------------------------------------------------
# MEANGS
# ---------------------------------------------------------------------------
def rename_scaffold_seeds(meangs_task_dir: Path, base_name: str) -> List[Path]:
    """Rename MEANGS scaffold_seeds.fa* to <base>_mg.fasta, <base>_mg_2.fasta, ..."""
    base = safe_name(base_name)
    targets = {"scaffold_seeds.fas", "scaffold_seeds.fa", "scaffold_seeds.fasta"}
    found = sorted([p for p in meangs_task_dir.rglob("*") if p.is_file() and p.name in targets])

    renamed: List[Path] = []
    for i, src in enumerate(found, start=1):
        new_name = f"{base}_mg.fasta" if i == 1 else f"{base}_mg_{i}.fasta"
        dest = meangs_task_dir / new_name
        j = 1
        while dest.exists():
            dest = meangs_task_dir / f"{base}_mg_{i}_{j}.fasta"
            j += 1
        try:
            src.replace(dest)
        except Exception:
            shutil.copy2(src, dest)
            try:
                src.unlink()
            except Exception:
                pass
        renamed.append(dest)
    return renamed


def clean_meangs_subdir(task_dir: Path, outname: str) -> int:
    """Remove MEANGS output subdir, keeping logs and renamed seed FASTAs."""
    d = task_dir / outname
    if d.exists() and d.is_dir():
        try:
            shutil.rmtree(d, ignore_errors=True)
            return 1
        except Exception:
            pass
    return 0


def run_meangs_task(
    task_id: str,
    fq1: Path,
    fq2: Path,
    meangs_root: Path,
    args,
) -> Tuple[bool, str]:
    tid = safe_name(task_id)
    task_dir = meangs_root / tid
    task_dir.mkdir(parents=True, exist_ok=True)

    outname = safe_name(args.meangs_outname)

    base_args = [
        "--silence",
        "-1", str(fq1),
        "-2", str(fq2),
        "-o", outname,
        "-t", str(args.threads_meangs),
        "-i", str(args.insert_meangs),
        "--species_class", str(args.species_class_meangs),
        "-n", str(args.nsample_meangs),
    ]

    for flag, attr in [
        ("--seqscaf",    "seqscaf_meangs"),
        ("--keepMinLen", "keepMinLen_meangs"),
    ]:
        val = getattr(args, attr, None)
        if val is not None and str(val).strip():
            base_args += [flag, str(val)]

    for flag, attr in [
        ("--deepin",     "deepin_meangs"),
        ("--clip",       "clip_meangs"),
        ("--keepIntMed", "keepIntMed_meangs"),
        ("--skipassem",  "skipassem_meangs"),
        ("--skipqc",     "skipqc_meangs"),
        ("--skiphmm",    "skiphmm_meangs"),
        ("--skipextend", "skipextend_meangs"),
    ]:
        if getattr(args, attr, False):
            base_args.append(flag)

    base_args += extra_args_from_list(getattr(args, "meangs_extra", []))

    cmd = build_tool_cmd(
        tool_label="MEANGS",
        exe_name=MEANGS_EXE_NAME,
        base_args=base_args,
        conda_env=getattr(args, "conda_env_meangs", None),
        tool_path=getattr(args, "meangs_path", None),
        conda_exe=getattr(args, "conda_exe", "conda"),
    )

    logging.info(f"[MEANGS] {task_id}")
    logging.info(f"[MEANGS_CWD] {task_dir}")
    logging.info(f"[MEANGS_CMD] {' '.join(cmd)}")

    proc = subprocess.run(cmd, cwd=str(task_dir), stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, text=True)

    (task_dir / "cmd_meangs.txt").write_text(
        f"(cwd={task_dir}) " + " ".join(cmd) + "\n", encoding="utf-8")
    (task_dir / "stdout_meangs.log").write_text(proc.stdout or "", encoding="utf-8",
                                                 errors="ignore")
    (task_dir / "stderr_meangs.log").write_text(proc.stderr or "", encoding="utf-8",
                                                 errors="ignore")
    (task_dir / "returncode_meangs.txt").write_text(str(proc.returncode) + "\n",
                                                     encoding="utf-8")

    renamed = rename_scaffold_seeds(task_dir, tid)
    if renamed:
        logging.info(f"[MEANGS_SEED] {task_id}: " + ", ".join(p.name for p in renamed))

    if getattr(args, "clean_meangs", True):
        if clean_meangs_subdir(task_dir, outname):
            logging.info(f"[CLEAN_MEANGS] {task_id}: removed '{outname}' subdir")

    if proc.returncode != 0:
        logging.error(f"[MEANGS_FAILED] {task_id} rc={proc.returncode}")
        logging.error("---- stdout tail ----\n" + tail_text(proc.stdout))
        logging.error("---- stderr tail ----\n" + tail_text(proc.stderr))
        return False, f"meangs_returncode={proc.returncode}"

    return True, "ok"


# ---------------------------------------------------------------------------
# Seed selection for GetOrganelle MT
# ---------------------------------------------------------------------------
def fasta_total_len(fa: Path) -> int:
    total = 0
    try:
        with fa.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.startswith(">"):
                    total += len(line.strip())
    except Exception:
        return 0
    return total


def pick_seed_fasta(
    meangs_task_dir: Path,
    task_id: str,
    seed_glob: Optional[str] = None,
) -> Optional[Path]:
    """Pick the best seed FASTA from a MEANGS output directory."""
    tid = safe_name(task_id)
    mg = sorted(meangs_task_dir.glob(f"{tid}_mg*.fasta"))
    if mg:
        return max(mg, key=fasta_total_len)
    if seed_glob:
        cands = [p for p in meangs_task_dir.glob(seed_glob) if p.is_file()]
        if cands:
            return max(cands, key=fasta_total_len)
    cands = []
    for pat in ("**/*.fa", "**/*.fasta", "**/*.fna"):
        cands.extend(p for p in meangs_task_dir.glob(pat) if p.is_file())
    if not cands:
        return None
    return max(cands, key=fasta_total_len)


# ---------------------------------------------------------------------------
# GetOrganelle FROM READS
# ---------------------------------------------------------------------------
def clean_one_get_dir(sample_get_dir: Path) -> int:
    """Remove large intermediate files/dirs produced by GetOrganelle-from-reads."""
    removed = 0
    for fn in _GETORG_CLEAN_FILES:
        p = sample_get_dir / fn
        if p.exists() and p.is_file():
            try:
                p.unlink()
                removed += 1
            except Exception:
                pass
    for dn in _GETORG_CLEAN_DIRS:
        d = sample_get_dir / dn
        if d.exists() and d.is_dir():
            try:
                shutil.rmtree(d, ignore_errors=True)
                removed += 1
            except Exception:
                pass
    if removed:
        logging.info(f"[CLEAN_GET] {sample_get_dir.name}: removed {removed} item(s)")
    return removed


def list_result_fastas_getreads(result_dir: Path) -> List[Path]:
    """
    Find GetOrganelle result FASTAs in *result_dir*.

    Prefers files matching 'complete' / 'circular' / 'final' patterns;
    falls back to any .fasta / .fa / .fna at the top level.
    Excludes known input file names (seed.fasta, label.fasta, genes.fasta).
    """
    if not result_dir.exists():
        return []

    priority_globs = [
        "*complete*.fasta", "*complete*.fa",   "*complete*.fna",
        "*circular*.fasta", "*final*.fasta",   "*result*.fasta",
        "*assembly*.fasta", "*.path_sequence.fasta",
    ]
    fallback_globs = ["*.fasta", "*.fa", "*.fna"]
    bad_names = {"seed.fasta", "label.fasta", "genes.fasta"}

    cand: List[Path] = []
    for pat in priority_globs + fallback_globs:
        cand.extend(sorted(result_dir.glob(pat)))

    seen: set = set()
    uniq: List[Path] = []
    for p in cand:
        if not p.is_file() or p.name.lower() in bad_names:
            continue
        key = str(p.resolve())
        if key not in seen:
            seen.add(key)
            uniq.append(p)

    def _sort_key(p: Path) -> Tuple[int, str]:
        n = p.name.lower()
        score = 9
        for kw, sc in [
            ("complete", 0), ("circular", 1), ("final", 2),
            ("result", 3),   ("path_sequence", 4),
        ]:
            if kw in n:
                score = min(score, sc)
        return (score, n)

    return sorted(uniq, key=_sort_key)


def fasta_has_circular(fa: Path, max_headers: int = 800) -> bool:
    """Return True if any FASTA header in *fa* contains 'circular' (case-insensitive)."""
    try:
        with fa.open("r", encoding="utf-8", errors="ignore") as f:
            n = 0
            for line in f:
                if line.startswith(">"):
                    n += 1
                    if "circular" in line.lower():
                        return True
                    if n >= max_headers:
                        break
    except Exception:
        return False
    return False


def run_getreads_task(
    *,
    task_id: str,
    fq1: Path,
    fq2: Path,
    out_dir: Path,
    meangs_root: Optional[Path],
    args,
    locus: str,  # "mt" or "nr"
) -> Tuple[bool, str]:
    """
    Run get_organelle_from_reads.py for a single sample and locus.

    For locus='mt':
      • Optionally uses a MEANGS-derived seed (use_meangs_seed=True).
      • Optionally uses a user-supplied seed (mt_seed_fasta).
      • MEANGS seed takes priority over mt_seed_fasta.

    For locus='nr':
      • REQUIRES nr_seed_fasta and nr_genes_fasta (GetOrganelle -F anonym).

    Returns (ok: bool, reason: str).
    """
    if locus not in ("mt", "nr"):
        return False, f"invalid_locus={locus}"

    tid = safe_name(task_id)
    get_dir = out_dir / tid
    get_dir.mkdir(parents=True, exist_ok=True)

    # ── MT seed resolution ──────────────────────────────────────────────────
    mt_seed: Optional[Path] = None
    if locus == "mt":
        # Priority 1: MEANGS-derived seed
        if getattr(args, "use_meangs_seed", False) and meangs_root is not None:
            mdir = meangs_root / tid
            if mdir.exists():
                mt_seed = pick_seed_fasta(
                    mdir, tid, seed_glob=getattr(args, "seed_glob_meangs", None)
                )
                if mt_seed:
                    logging.info(f"[SEED_MT] {task_id} using MEANGS seed: {mt_seed}")
                else:
                    logging.info(f"[SEED_MT] {task_id}: no MEANGS seed found")

        # Priority 2: user-supplied mt_seed_fasta
        if mt_seed is None:
            user_seed = getattr(args, "mt_seed_fasta", None)
            if user_seed:
                p = Path(str(user_seed))
                if p.is_file():
                    mt_seed = p
                    logging.info(f"[SEED_MT] {task_id} using user seed: {mt_seed}")
                else:
                    logging.warning(
                        f"[SEED_MT] {task_id}: --mt_seed_fasta not found: {p} (skipping -s)"
                    )

    # ── NR seed / genes validation ──────────────────────────────────────────
    if locus == "nr":
        nr_seed  = getattr(args, "nr_seed_fasta",  None)
        nr_genes = getattr(args, "nr_genes_fasta", None)
        if not nr_seed or not Path(str(nr_seed)).is_file():
            return False, "nr_seed_fasta_missing_or_not_a_file"
        if not nr_genes or not Path(str(nr_genes)).is_file():
            return False, "nr_genes_fasta_missing_or_not_a_file"

    # ── Build argument list ─────────────────────────────────────────────────
    if locus == "mt":
        base_args = [
            "-1", str(fq1),
            "-2", str(fq2),
            "-o", str(get_dir),
            "-F", str(getattr(args, "F_mt", "animal_mt")),
            "-R", str(getattr(args, "R_mt",  10)),
            "-t", str(getattr(args, "t_mt",   8)),
            "-k", str(getattr(args, "k_mt",  "21,45,65,85,105")),
            "--continue",
        ]
        if mt_seed:
            base_args += ["-s", str(mt_seed)]

        # Optional word size
        w_mt = getattr(args, "w_mt", None)
        if w_mt is not None and str(w_mt).strip():
            base_args += ["-w", str(w_mt)]

        # Optional read cap (speeds up large datasets)
        max_reads_mt = getattr(args, "max_reads_mt", None)
        if max_reads_mt is not None and str(max_reads_mt).strip():
            base_args += ["--max-reads", str(max_reads_mt)]

        if getattr(args, "all_data_mt", False):
            base_args += ["--max-reads", "inf", "--reduce-reads-for-coverage", "inf"]

        base_args += extra_args_from_list(getattr(args, "get_extra_mt", []))

    else:  # nr
        base_args = [
            "-1", str(fq1),
            "-2", str(fq2),
            "-o", str(get_dir),
            "-F", str(getattr(args, "F_nr", "anonym")),
            "-R", str(getattr(args, "R_nr",  10)),
            "-t", str(getattr(args, "t_nr",   8)),
            "-k", str(getattr(args, "k_nr",  "21,45,65,85,105")),
            "--continue",
            "-s",      str(getattr(args, "nr_seed_fasta")),
            "--genes", str(getattr(args, "nr_genes_fasta")),
        ]
        p_nr = getattr(args, "P_nr", None)
        if p_nr is not None and str(p_nr).strip():
            base_args += ["-P", str(p_nr)]

        maxext = getattr(args, "max_extending_len_nr", None)
        if maxext is not None and str(maxext).strip():
            base_args += ["--max-extending-len", str(maxext)]

        max_reads_nr = getattr(args, "max_reads_nr", None)
        if max_reads_nr is not None and str(max_reads_nr).strip():
            base_args += ["--max-reads", str(max_reads_nr)]

        if getattr(args, "all_data_nr", False):
            base_args += ["--max-reads", "inf", "--reduce-reads-for-coverage", "inf"]

        base_args += extra_args_from_list(getattr(args, "get_extra_nr", []))

    cmd = build_tool_cmd(
        tool_label="GetOrganelle",
        exe_name=GETORG_READS_EXE_NAME,
        base_args=base_args,
        conda_env=getattr(args, "conda_env_getorganelle", None),
        tool_path=getattr(args, "getorganelle_path", None),
        conda_exe=getattr(args, "conda_exe", "conda"),
    )

    logging.info(f"[GET_READS] {task_id} locus={locus}")
    logging.info(f"[GET_CMD]   {' '.join(cmd)}")

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    (get_dir / f"cmd_getorganelle_{locus}.txt").write_text(
        " ".join(cmd) + "\n", encoding="utf-8")
    (get_dir / f"stdout_getorganelle_{locus}.log").write_text(
        proc.stdout or "", encoding="utf-8", errors="ignore")
    (get_dir / f"stderr_getorganelle_{locus}.log").write_text(
        proc.stderr or "", encoding="utf-8", errors="ignore")
    (get_dir / f"returncode_getorganelle_{locus}.txt").write_text(
        str(proc.returncode) + "\n", encoding="utf-8")

    if proc.returncode != 0:
        logging.error(f"[GET_FAILED] {task_id} locus={locus} rc={proc.returncode}")
        logging.error("---- stdout tail ----\n" + tail_text(proc.stdout))
        logging.error("---- stderr tail ----\n" + tail_text(proc.stderr))
        return False, f"getorganelle_returncode={proc.returncode}"

    fastas = list_result_fastas_getreads(get_dir)
    if not fastas:
        logging.error(f"[GET_NO_FASTA] {task_id} locus={locus}: assembly completed but no FASTA found")
        return False, "no_fasta_output"

    if getattr(args, "clean_get", False):
        clean_one_get_dir(get_dir)

    return True, "ok"


# ---------------------------------------------------------------------------
# Summary & reporting
# ---------------------------------------------------------------------------
def write_renamed_fasta(
    src_fa: Path,
    out_fh,
    sample: str,
    locus: str,
    file_idx: int,
    topology: str,
) -> Tuple[int, int, int]:
    """
    Write normalised FASTA to *out_fh* with standardised headers:
      >{sample}|locus={locus}|{file_idx} topology={topology}

    Returns (n_seqs, total_bp, n_N).
    """
    seq_idx = 0
    nseq = total = nN = 0
    last = ""
    with src_fa.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            last = line
            if line.startswith(">"):
                seq_idx += 1
                nseq += 1
                tag = (
                    f"{sample}|locus={locus}|{file_idx}"
                    if seq_idx == 1
                    else f"{sample}|locus={locus}|{file_idx}_{seq_idx}"
                )
                out_fh.write(f">{tag} topology={topology}\n")
            else:
                s = line.strip()
                if s:
                    total += len(s)
                    nN += s.upper().count("N")
                out_fh.write(line)
    if last and not last.endswith("\n"):
        out_fh.write("\n")
    out_fh.write("\n")
    return nseq, total, nN


def check_one_locus(
    out_root: Path,
    dir_name: str,
    expected_task_ids: List[str],
) -> List[str]:
    droot = out_root / dir_name
    failed: List[str] = []
    for tid in expected_task_ids:
        sdir = droot / tid
        if not sdir.exists() or not sdir.is_dir() or not list_result_fastas_getreads(sdir):
            failed.append(tid)
    return failed


def write_list(names: List[str], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    uniq = sorted(set(names))
    out_path.write_text("\n".join(uniq) + ("\n" if uniq else ""), encoding="utf-8")
    logging.info(f"[WRITE] {out_path} (n={len(uniq)})")


def summarize_one_locus(
    out_root: Path,
    dir_name: str,
    locus: str,
    fasta_out: Path,
    report_rows: List[List[str]],
) -> None:
    droot = out_root / dir_name
    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    with fasta_out.open("w", encoding="utf-8") as out:
        for sdir in sorted([p for p in droot.iterdir() if p.is_dir()]):
            sample = sdir.name
            fastas = list_result_fastas_getreads(sdir)
            if not fastas:
                report_rows.append([sample, locus, "0", "", "", "0", "0", "0", str(sdir)])
                continue
            for i, fa in enumerate(fastas, start=1):
                topo = "circular" if fasta_has_circular(fa) else "linear"
                nseq, total, nN = write_renamed_fasta(
                    fa, out, sample=sample, locus=locus, file_idx=i, topology=topo
                )
                report_rows.append(
                    [sample, locus, "1", str(fa), topo, str(nseq), str(total), str(nN), str(sdir)]
                )


def write_report_tsv(report_path: Path, rows: List[List[str]]) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "sample", "locus", "has_fasta", "fasta_path",
        "topology", "n_seqs", "total_bp", "n_N", "sample_outdir",
    ]
    with report_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    logging.info(f"[REPORT] {report_path} rows={len(rows)}")


# ---------------------------------------------------------------------------
# Memory detection
# ---------------------------------------------------------------------------
def detect_total_mem_mb_linux() -> Optional[int]:
    p = Path("/proc/meminfo")
    if not p.exists():
        return None
    try:
        for line in p.read_text(encoding="utf-8", errors="ignore").splitlines():
            if line.startswith("MemTotal:"):
                parts = line.split()
                if len(parts) >= 2 and parts[1].isdigit():
                    return int(parts[1]) // 1024
    except Exception:
        return None
    return None


# ---------------------------------------------------------------------------
# Snakefile templates
# ---------------------------------------------------------------------------
_SNAKEFILE_HEADER = '''configfile: "config.json"
import os

DRIVER  = config["driver_py"]
PY      = config["python_exe"]
SAMPLES = sorted(list(config.get("samples", {}).keys()))
'''

SNAKEFILE_MGS2GET = (
    _SNAKEFILE_HEADER
    + r"""
MEANGS_DONE = os.path.join(config["out_root"], config["meangs_dir"], "{sample}", ".meangs.done")
MT_DONE     = os.path.join(config["out_root"], config["mt_dir"],     "{sample}", ".mt.done")

rule all:
    input:
        config["summary_all"], config["failed_all"], config["report_tsv"]

rule meangs:
    input:
        r1=lambda wc: config["samples"][wc.sample]["r1"],
        r2=lambda wc: config["samples"][wc.sample]["r2"]
    output: MEANGS_DONE
    threads: int(config["params"]["threads_meangs"])
    resources:
        mem_mb=int(config["params"]["mem_meangs_mb"]),
        meangs_slots=1
    shell:
        f'"{PY}" "{DRIVER}" --smk_internal meangs_one --config config.json --sample {{wildcards.sample}}'

rule mt:
    input:
        r1=lambda wc: config["samples"][wc.sample]["r1"],
        r2=lambda wc: config["samples"][wc.sample]["r2"],
        m_done=MEANGS_DONE
    output: MT_DONE
    threads: int(config["params"]["t_mt"])
    resources:
        mem_mb=int(config["params"]["mem_mt_mb"]),
        mt_slots=1
    shell:
        f'"{PY}" "{DRIVER}" --smk_internal mt_one --config config.json --sample {{wildcards.sample}}'

rule finalize:
    input:  expand(MT_DONE, sample=SAMPLES)
    output:
        config["report_tsv"], config["failed_mt"],
        config["failed_all"], config["summary_mt"], config["summary_all"]
    shell:
        f'"{PY}" "{DRIVER}" --smk_internal finalize_mt --config config.json'
"""
)

SNAKEFILE_MT_ONLY = (
    _SNAKEFILE_HEADER
    + r"""
MT_DONE = os.path.join(config["out_root"], config["mt_dir"], "{sample}", ".mt.done")

rule all:
    input:
        config["summary_all"], config["failed_all"], config["report_tsv"]

rule mt:
    input:
        r1=lambda wc: config["samples"][wc.sample]["r1"],
        r2=lambda wc: config["samples"][wc.sample]["r2"]
    output: MT_DONE
    threads: int(config["params"]["t_mt"])
    resources:
        mem_mb=int(config["params"]["mem_mt_mb"]),
        mt_slots=1
    shell:
        f'"{PY}" "{DRIVER}" --smk_internal mt_one --config config.json --sample {{wildcards.sample}}'

rule finalize:
    input:  expand(MT_DONE, sample=SAMPLES)
    output:
        config["report_tsv"], config["failed_mt"],
        config["failed_all"], config["summary_mt"], config["summary_all"]
    shell:
        f'"{PY}" "{DRIVER}" --smk_internal finalize_mt --config config.json'
"""
)

SNAKEFILE_NR_ONLY = (
    _SNAKEFILE_HEADER
    + r"""
NR_DONE = os.path.join(config["out_root"], config["nr_dir"], "{sample}", ".nr.done")

rule all:
    input:
        config["summary_all"], config["failed_all"], config["report_tsv"]

rule nr:
    input:
        r1=lambda wc: config["samples"][wc.sample]["r1"],
        r2=lambda wc: config["samples"][wc.sample]["r2"]
    output: NR_DONE
    threads: int(config["params"]["t_nr"])
    resources:
        mem_mb=int(config["params"]["mem_nr_mb"]),
        nr_slots=1
    shell:
        f'"{PY}" "{DRIVER}" --smk_internal nr_one --config config.json --sample {{wildcards.sample}}'

rule finalize:
    input:  expand(NR_DONE, sample=SAMPLES)
    output:
        config["report_tsv"], config["failed_nr"],
        config["failed_all"], config["summary_nr"], config["summary_all"]
    shell:
        f'"{PY}" "{DRIVER}" --smk_internal finalize_nr --config config.json'
"""
)

SNAKEFILE_MT_NR = (
    _SNAKEFILE_HEADER
    + r"""
MT_DONE = os.path.join(config["out_root"], config["mt_dir"], "{sample}", ".mt.done")
NR_DONE = os.path.join(config["out_root"], config["nr_dir"], "{sample}", ".nr.done")

rule all:
    input:
        config["summary_all"], config["failed_all"], config["report_tsv"]

rule mt:
    input:
        r1=lambda wc: config["samples"][wc.sample]["r1"],
        r2=lambda wc: config["samples"][wc.sample]["r2"]
    output: MT_DONE
    threads: int(config["params"]["t_mt"])
    resources:
        mem_mb=int(config["params"]["mem_mt_mb"]),
        mt_slots=1
    shell:
        f'"{PY}" "{DRIVER}" --smk_internal mt_one --config config.json --sample {{wildcards.sample}}'

rule nr:
    input:
        r1=lambda wc: config["samples"][wc.sample]["r1"],
        r2=lambda wc: config["samples"][wc.sample]["r2"]
    output: NR_DONE
    threads: int(config["params"]["t_nr"])
    resources:
        mem_mb=int(config["params"]["mem_nr_mb"]),
        nr_slots=1
    shell:
        f'"{PY}" "{DRIVER}" --smk_internal nr_one --config config.json --sample {{wildcards.sample}}'

rule finalize:
    input:
        expand(MT_DONE, sample=SAMPLES),
        expand(NR_DONE, sample=SAMPLES)
    output:
        config["report_tsv"], config["failed_mt"], config["failed_nr"],
        config["failed_all"], config["summary_mt"], config["summary_nr"], config["summary_all"]
    shell:
        f'"{PY}" "{DRIVER}" --smk_internal finalize_mt_nr --config config.json'
"""
)

SNAKEFILE_ALL = (
    _SNAKEFILE_HEADER
    + r"""
MEANGS_DONE = os.path.join(config["out_root"], config["meangs_dir"], "{sample}", ".meangs.done")
MT_DONE     = os.path.join(config["out_root"], config["mt_dir"],     "{sample}", ".mt.done")
NR_DONE     = os.path.join(config["out_root"], config["nr_dir"],     "{sample}", ".nr.done")

rule all:
    input:
        config["report_tsv"], config["failed_mt"], config["failed_nr"],
        config["failed_all"], config["summary_mt"], config["summary_nr"], config["summary_all"]

rule meangs:
    input:
        r1=lambda wc: config["samples"][wc.sample]["r1"],
        r2=lambda wc: config["samples"][wc.sample]["r2"]
    output: MEANGS_DONE
    threads: int(config["params"]["threads_meangs"])
    resources:
        mem_mb=int(config["params"]["mem_meangs_mb"]),
        meangs_slots=1
    shell:
        f'"{PY}" "{DRIVER}" --smk_internal meangs_one --config config.json --sample {{wildcards.sample}}'

rule mt:
    input:
        r1=lambda wc: config["samples"][wc.sample]["r1"],
        r2=lambda wc: config["samples"][wc.sample]["r2"],
        m_done=MEANGS_DONE
    output: MT_DONE
    threads: int(config["params"]["t_mt"])
    resources:
        mem_mb=int(config["params"]["mem_mt_mb"]),
        mt_slots=1
    shell:
        f'"{PY}" "{DRIVER}" --smk_internal mt_one --config config.json --sample {{wildcards.sample}}'

rule nr:
    input:
        r1=lambda wc: config["samples"][wc.sample]["r1"],
        r2=lambda wc: config["samples"][wc.sample]["r2"],
        mt_done=MT_DONE
    output: NR_DONE
    threads: int(config["params"]["t_nr"])
    resources:
        mem_mb=int(config["params"]["mem_nr_mb"]),
        nr_slots=1
    shell:
        f'"{PY}" "{DRIVER}" --smk_internal nr_one --config config.json --sample {{wildcards.sample}}'

rule finalize:
    input:
        expand(MT_DONE, sample=SAMPLES),
        expand(NR_DONE, sample=SAMPLES)
    output:
        config["report_tsv"], config["failed_mt"], config["failed_nr"],
        config["failed_all"], config["summary_mt"], config["summary_nr"], config["summary_all"]
    shell:
        f'"{PY}" "{DRIVER}" --smk_internal finalize_mt_nr --config config.json'
"""
)


def write_smk_files(workflow_dir: Path, config: dict, mode: str) -> Tuple[Path, Path]:
    workflow_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = workflow_dir / "config.json"
    sf_path  = workflow_dir / "Snakefile"
    cfg_path.write_text(json.dumps(config, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    _map = {
        "mgs2get":    SNAKEFILE_MGS2GET,
        "mt-only":    SNAKEFILE_MT_ONLY,
        "nr-only":    SNAKEFILE_NR_ONLY,
        "mt-nr":      SNAKEFILE_MT_NR,
        "all":        SNAKEFILE_ALL,
    }
    if mode not in _map:
        raise ValueError(f"Unknown mode: {mode!r}.  Valid: {list(_map)}")
    sf_path.write_text(_map[mode].lstrip() + "\n", encoding="utf-8")
    return sf_path, cfg_path


def run_snakemake_cli(
    *,
    workflow_dir: Path,
    snakefile: Path,
    cores: int,
    jobs: int,
    latency_wait: int,
    dryrun: bool,
    mem_total_mb: Optional[int],
    max_tasks_meangs: int,
    max_tasks_mt: int,
    max_tasks_nr: int,
) -> int:
    smk = which_or_none("snakemake")
    if not smk:
        raise FileNotFoundError(
            "snakemake not found in PATH.\n"
            "Install via:  conda install -c bioconda snakemake"
        )
    cmd = [
        smk,
        "-s", str(snakefile),
        "--cores",           str(cores),
        "--jobs",            str(jobs),
        "--rerun-incomplete",
        "--keep-going",
        "--printshellcmds",
        "--latency-wait",    str(latency_wait),
    ]

    res_parts: List[str] = []
    if mem_total_mb is not None:
        res_parts.append(f"mem_mb={int(mem_total_mb)}")
    res_parts += [
        f"meangs_slots={int(max_tasks_meangs)}",
        f"mt_slots={int(max_tasks_mt)}",
        f"nr_slots={int(max_tasks_nr)}",
    ]
    cmd += ["--resources"] + res_parts

    if dryrun:
        cmd += ["-n"]

    logging.info("[SNAKEMAKE_CMD] " + " ".join(cmd))
    p = subprocess.run(cmd, cwd=str(workflow_dir))
    return int(p.returncode)


# ---------------------------------------------------------------------------
# Internal Snakemake workers  (called via --smk_internal)
# ---------------------------------------------------------------------------
def load_config(config_path: Path) -> dict:
    return json.loads(config_path.read_text(encoding="utf-8"))


def ns_from_params(params: dict) -> SimpleNamespace:
    return SimpleNamespace(**params)


def _stamp_write(stamp: Path, ok: bool, reason: str) -> None:
    stamp.parent.mkdir(parents=True, exist_ok=True)
    stamp.write_text(f"ok={1 if ok else 0}\nreason={reason}\n", encoding="utf-8")


def _smk_meangs_one(config_path: Path, sample: str) -> int:
    cfg    = load_config(config_path)
    params = ns_from_params(cfg["params"])
    sample = safe_name(sample)

    out_root    = Path(cfg["out_root"]).resolve()
    meangs_root = out_root / cfg["meangs_dir"]
    task_dir    = meangs_root / sample
    task_dir.mkdir(parents=True, exist_ok=True)
    setup_logger(task_dir / "task_meangs.log")

    s = cfg["samples"].get(sample)
    if not s:
        _stamp_write(task_dir / ".meangs.done", False, "sample_not_in_config")
        return 0

    fq1 = Path(s["r1"]).resolve()
    fq2 = Path(s["r2"]).resolve()
    ok, reason = False, "unknown"
    try:
        ok, reason = run_meangs_task(sample, fq1, fq2, meangs_root, params)
    except Exception as e:
        logging.error(f"[MEANGS_EXCEPTION] {sample}: {e}")
        ok, reason = False, f"exception:{e}"

    _stamp_write(task_dir / ".meangs.done", ok, reason)
    return 0


def _smk_mt_one(config_path: Path, sample: str) -> int:
    cfg    = load_config(config_path)
    params = ns_from_params(cfg["params"])
    sample = safe_name(sample)

    out_root    = Path(cfg["out_root"]).resolve()
    mt_root     = out_root / cfg["mt_dir"]
    meangs_root = out_root / cfg["meangs_dir"]
    mt_dir      = mt_root / sample
    mt_dir.mkdir(parents=True, exist_ok=True)
    setup_logger(mt_dir / "task_mt.log")

    s = cfg["samples"].get(sample)
    if not s:
        _stamp_write(mt_dir / ".mt.done", False, "sample_not_in_config")
        return 0

    fq1 = Path(s["r1"]).resolve()
    fq2 = Path(s["r2"]).resolve()
    ok, reason = False, "unknown"
    try:
        ok, reason = run_getreads_task(
            task_id=sample,
            fq1=fq1,
            fq2=fq2,
            out_dir=mt_root,
            meangs_root=meangs_root if cfg.get("use_meangs_for_mt", False) else None,
            args=params,
            locus="mt",
        )
    except Exception as e:
        logging.error(f"[MT_EXCEPTION] {sample}: {e}")
        ok, reason = False, f"exception:{e}"

    _stamp_write(mt_dir / ".mt.done", ok, reason)
    return 0


def _smk_nr_one(config_path: Path, sample: str) -> int:
    cfg    = load_config(config_path)
    params = ns_from_params(cfg["params"])
    sample = safe_name(sample)

    out_root = Path(cfg["out_root"]).resolve()
    nr_root  = out_root / cfg["nr_dir"]
    nr_dir   = nr_root / sample
    nr_dir.mkdir(parents=True, exist_ok=True)
    setup_logger(nr_dir / "task_nr.log")

    s = cfg["samples"].get(sample)
    if not s:
        _stamp_write(nr_dir / ".nr.done", False, "sample_not_in_config")
        return 0

    fq1 = Path(s["r1"]).resolve()
    fq2 = Path(s["r2"]).resolve()
    ok, reason = False, "unknown"
    try:
        ok, reason = run_getreads_task(
            task_id=sample,
            fq1=fq1,
            fq2=fq2,
            out_dir=nr_root,
            meangs_root=None,
            args=params,
            locus="nr",
        )
    except Exception as e:
        logging.error(f"[NR_EXCEPTION] {sample}: {e}")
        ok, reason = False, f"exception:{e}"

    _stamp_write(nr_dir / ".nr.done", ok, reason)
    return 0


def _finalize_common(cfg: dict, which: str) -> int:
    out_root = Path(cfg["out_root"]).resolve()
    setup_logger(out_root / "pipeline.log")

    samples       = sorted(list(cfg.get("samples", {}).keys()))
    missing_reads = list(cfg.get("missing_reads", []))
    rows:     List[List[str]] = []
    mt_failed: List[str] = []
    nr_failed: List[str] = []

    if which in ("mt", "mt-nr"):
        mt_failed = check_one_locus(out_root, cfg["mt_dir"], samples)
        write_list(sorted(set(mt_failed) | set(missing_reads)), Path(cfg["failed_mt"]))
        summarize_one_locus(out_root, cfg["mt_dir"], "mt", Path(cfg["summary_mt"]), rows)

    if which in ("nr", "mt-nr"):
        nr_failed = check_one_locus(out_root, cfg["nr_dir"], samples)
        write_list(sorted(set(nr_failed) | set(missing_reads)), Path(cfg["failed_nr"]))
        summarize_one_locus(out_root, cfg["nr_dir"], "nr", Path(cfg["summary_nr"]), rows)

    failed_all = set(missing_reads)
    if which in ("mt", "mt-nr"):
        failed_all |= set(mt_failed)
    if which in ("nr", "mt-nr"):
        failed_all |= set(nr_failed)
    write_list(sorted(failed_all), Path(cfg["failed_all"]))

    summary_all = Path(cfg["summary_all"])
    summary_all.parent.mkdir(parents=True, exist_ok=True)
    with summary_all.open("w", encoding="utf-8") as out:
        for key in (("mt", "summary_mt"), ("nr", "summary_nr")):
            locus, k = key
            if (which in (locus, "mt-nr")) and k in cfg:
                p = Path(cfg[k])
                if p.exists():
                    out.write(p.read_text(encoding="utf-8", errors="ignore"))

    write_report_tsv(Path(cfg["report_tsv"]), rows)
    return 0


def _smk_finalize_mt(config_path: Path)    -> int: return _finalize_common(load_config(config_path), "mt")
def _smk_finalize_nr(config_path: Path)    -> int: return _finalize_common(load_config(config_path), "nr")
def _smk_finalize_mt_nr(config_path: Path) -> int: return _finalize_common(load_config(config_path), "mt-nr")


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------
class _HelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Slightly wider help formatter."""
    def __init__(self, prog, **kw):
        super().__init__(prog, max_help_position=36, width=100, **kw)


def _add_common_io(sp: argparse.ArgumentParser) -> None:
    g = sp.add_argument_group("I/O paths")
    g.add_argument("-o", "--out_root",   default="mgs2get_results",
                   help="Output root directory  (default: %(default)s)")
    g.add_argument("--meangs_dir",       default="01.MEANGS_results",
                   help="MEANGS results subdir  (default: %(default)s)")
    g.add_argument("--mt_dir",           default="02_mt",
                   help="MT GetOrganelle results subdir  (default: %(default)s)")
    g.add_argument("--nr_dir",           default="03_nr",
                   help="NR/rRNA GetOrganelle results subdir  (default: %(default)s)")


def _add_reads_layout(sp: argparse.ArgumentParser) -> None:
    g = sp.add_argument_group("Reads layout / scanning")
    g.add_argument("-i", "--input_dir", required=True,
                   help="Directory containing raw reads (subdirs or flat layout).")
    g.add_argument(
        "--fq_position",
        choices=["subdir", "flat"],
        default="subdir",
        metavar="subdir|flat",
        help=textwrap.dedent("""\
            How reads are organised:
              subdir  input_dir/SAMPLE/*_1.fq.gz  (default)
              flat    input_dir/*_1.fq.gz"""),
    )
    g.add_argument(
        "--suffix_fq",
        nargs="+",
        default=DEFAULT_SUFFIX_PAIRS,
        metavar="R1,R2[;R1,R2...]",
        help=textwrap.dedent("""\
            Paired-end read suffix pairs, tried in order (first match wins).
            Formats:
              single pair        "_1.fq.gz,_2.fq.gz"
              semicolon-joined   "_1.clean.fq.gz,_2.clean.fq.gz;_1.fq.gz,_2.fq.gz"
              space-separated    "_1.clean.fq.gz,_2.clean.fq.gz" "_1.fq.gz,_2.fq.gz"
            Default tries 5 common patterns automatically."""),
    )
    g.add_argument(
        "--subdir_list",
        nargs="*",
        default=None,
        metavar="SAMPLE",
        help="Whitelist of sample names / subdir names.  "
             "Accepts: space-separated values, comma-separated, or a .txt file (one per line).",
    )


def _add_exec_control(sp: argparse.ArgumentParser) -> None:
    g = sp.add_argument_group("Tool execution")
    g.add_argument("--conda_exe",               default="conda",
                   help="conda/mamba executable  (default: %(default)s)")
    g.add_argument("--conda_env_meangs",         default="meangs",
                   help="Conda env for MEANGS  (default: %(default)s)")
    g.add_argument("--conda_env_getorganelle",   default="getorganelle",
                   help="Conda env for GetOrganelle  (default: %(default)s)")
    g.add_argument("--meangs_path",              default=None,
                   help="Path to meangs.py or its parent directory")
    g.add_argument("--getorganelle_path",        default=None,
                   help="Path to get_organelle_from_reads.py or its parent directory")


def _add_meangs_params(sp: argparse.ArgumentParser) -> None:
    g = sp.add_argument_group("MEANGS parameters")
    g.add_argument("--meangs_outname",       default="meangs",
                   help="MEANGS -o value (relative dir, no slashes)  (default: %(default)s)")
    g.add_argument("--threads_meangs",       type=int, default=4,    metavar="N")
    g.add_argument("--insert_meangs",        type=int, default=350,  metavar="N",
                   help="Insert size estimate  (default: %(default)s)")
    g.add_argument("--species_class_meangs", default="A-worms",
                   help="MEANGS --species_class  (default: %(default)s)")
    g.add_argument("--nsample_meangs",       type=int, default=2000000, metavar="N")
    g.add_argument("--seqscaf_meangs",       default=None,   metavar="FILE")
    g.add_argument("--keepMinLen_meangs",    type=int, default=None,  metavar="N")
    g.add_argument("--deepin_meangs",        action="store_true", default=True,
                   help="MEANGS --deepin  (default ON)")
    g.add_argument("--no_deepin_meangs",     action="store_false", dest="deepin_meangs",
                   help="Disable --deepin")
    g.add_argument("--clip_meangs",          action="store_true")
    g.add_argument("--keepIntMed_meangs",    action="store_true")
    g.add_argument("--skipassem_meangs",     action="store_true")
    g.add_argument("--skipqc_meangs",        action="store_true")
    g.add_argument("--skiphmm_meangs",       action="store_true")
    g.add_argument("--skipextend_meangs",    action="store_true")
    g.add_argument("--meangs_extra",         action="append", default=[], metavar="ARG",
                   help="Extra arguments passed verbatim to MEANGS (repeatable)")


def _add_mt_params(sp: argparse.ArgumentParser) -> None:
    g = sp.add_argument_group("GetOrganelle MT assembly parameters")
    g.add_argument("--F_mt",  default="animal_mt",
                   help="GetOrganelle -F for MT  (default: %(default)s)\n"
                        "Other options: embplant_pt, embplant_mt, fungus_mt, ...")
    g.add_argument("--R_mt",  type=int, default=10,            metavar="N",
                   help="Max extension rounds  (default: %(default)s; fewer = faster)")
    g.add_argument("--t_mt",  type=int, default=8,             metavar="N",
                   help="Threads per MT job  (default: %(default)s)")
    g.add_argument("--k_mt",  default="21,45,65,85,105",
                   help="k-mer sizes for SPAdes  (default: %(default)s)")
    g.add_argument(
        "--mt_seed_fasta",
        default=None,
        metavar="FILE",
        help=textwrap.dedent("""\
            OPTIONAL: Provide a seed sequence for MT assembly (-s).
            Using a seed from a closely related species dramatically SPEEDS UP
            GetOrganelle by guiding the initial read-filtering step.
            Examples:
              --mt_seed_fasta /db/Mytilus_edulis_mt.fasta
              --mt_seed_fasta /db/animal_mt.fasta
            If --use_meangs_seed is also enabled, the MEANGS-derived seed takes
            priority and this value is used as a fallback."""),
    )
    g.add_argument("--w_mt",  default=None, metavar="INT",
                   help="Word size for GetOrganelle graph (-w).  Larger = faster but less "
                        "sensitive.  Leave blank for auto-detection (recommended).")
    g.add_argument("--max_reads_mt", default=None, metavar="INT",
                   help="Cap reads used for MT assembly (--max-reads).  "
                        "1500000 is often sufficient.  Default: no cap.")
    g.add_argument("--all_data_mt",  action="store_true",
                   help="Use ALL reads (overrides --max_reads_mt; sets --max-reads inf)")

    g2 = sp.add_argument_group("MT seed from MEANGS output (mgs2get / all modes)")
    g2.add_argument("--use_meangs_seed",    action="store_true",  default=True,
                    help="Use MEANGS scaffold_seeds as GetOrganelle -s seed  (default ON)")
    g2.add_argument("--no_use_meangs_seed", action="store_false", dest="use_meangs_seed",
                    help="Disable MEANGS seed")
    g2.add_argument("--seed_glob_meangs",   default=None, metavar="GLOB",
                    help="Glob pattern to find seed FASTA in MEANGS output dir "
                         "(e.g. '*_mg.fasta')")

    g.add_argument("--get_extra_mt", action="append", default=[], metavar="ARG",
                   help="Extra arguments passed verbatim to GetOrganelle MT (repeatable)")


def _add_nr_params(sp: argparse.ArgumentParser) -> None:
    g = sp.add_argument_group(
        "GetOrganelle NR/rRNA assembly  (-F anonym + -s + --genes)",
        description=textwrap.dedent("""\
            GetOrganelle assembles custom loci (18S, 28S, ITS, etc.) using
            '-F anonym' which REQUIRES both a seed FASTA (-s) and a gene/label
            FASTA (--genes).  See GetOrganelle wiki:
            https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#custom-loci"""),
    )
    g.add_argument("--F_nr",  default="anonym",
                   help="GetOrganelle -F for NR  (default: %(default)s; do not change)")
    g.add_argument(
        "--nr_seed_fasta",
        default=None,
        metavar="FILE",
        required=False,
        help="REQUIRED for NR modes: seed FASTA for -s.\n"
             "Typically concatenated 18S/28S reference sequences from NCBI.\n"
             "Example: /db/rRNA/18-28s_seed.fasta",
    )
    g.add_argument(
        "--nr_genes_fasta",
        default=None,
        metavar="FILE",
        required=False,
        help="REQUIRED for NR modes: label/gene FASTA for --genes.\n"
             "Used to guide graph disentanglement.\n"
             "Example: /db/rRNA/18-28s_label.fasta",
    )
    g.add_argument("--R_nr",  type=int, default=10,           metavar="N")
    g.add_argument("--t_nr",  type=int, default=8,            metavar="N")
    g.add_argument("--k_nr",  default="21,45,65,85,105",
                   help="k-mer sizes for NR SPAdes  (default: %(default)s)")
    g.add_argument("--P_nr",  default="0",
                   help="GetOrganelle -P (path score cutoff).  "
                        "0 = keep all paths; RECOMMENDED for 18S/28S.  (default: %(default)s)")
    g.add_argument("--max_extending_len_nr", default="",  metavar="INT",
                   help="GetOrganelle --max-extending-len.  "
                        "250 is a good starting point for rRNA (~1800 bp).  "
                        "Leave blank for software default.")
    g.add_argument("--max_reads_nr", default=None, metavar="INT",
                   help="Cap reads for NR assembly (--max-reads).  Default: no cap.")
    g.add_argument("--all_data_nr",  action="store_true",
                   help="Use ALL reads for NR  (sets --max-reads inf)")
    g.add_argument("--get_extra_nr", action="append", default=[], metavar="ARG",
                   help="Extra arguments passed verbatim to GetOrganelle NR (repeatable)")


def _add_outputs(sp: argparse.ArgumentParser) -> None:
    g = sp.add_argument_group("Output file names & cleanup")
    g.add_argument("--clean_meangs",    action="store_true",  default=True,
                   help="Remove MEANGS raw output subdir (keep logs+seeds)  (default ON)")
    g.add_argument("--no_clean_meangs", action="store_false", dest="clean_meangs")
    g.add_argument("--clean_get",       action="store_true",  default=True,
                   help="Remove GetOrganelle large temp dirs  (default ON)")
    g.add_argument("--no_clean_get",    action="store_false", dest="clean_get")

    g.add_argument("--summary_mt_name",  default="summary_fasta_mt.txt")
    g.add_argument("--summary_nr_name",  default="summary_fasta_nr.txt")
    g.add_argument("--summary_all_name", default="summary_fasta_all.txt")
    g.add_argument("--report_name",      default="summary_report.tsv")
    g.add_argument("--failed_mt_name",   default="mt_failed.txt")
    g.add_argument("--failed_nr_name",   default="nr_failed.txt")
    g.add_argument("--failed_all_name",  default="failed_all.txt")


def _add_smk_sched(sp: argparse.ArgumentParser) -> None:
    g = sp.add_argument_group("Snakemake scheduling")
    g.add_argument("--max_parallel",     type=int, default=8,
                   help="Max concurrent sample jobs  (snakemake --jobs)  (default: %(default)s)")
    g.add_argument("--smk_cores",        type=int, default=0,
                   help="Total CPU budget  (snakemake --cores).  0 = auto  (default: %(default)s)")
    g.add_argument("--latency_wait",     type=int, default=60,
                   help="Snakemake --latency-wait (s)  (default: %(default)s)")
    g.add_argument("--dryrun",           action="store_true",
                   help="Show workflow plan without running  (snakemake -n)")
    g.add_argument("--max_tasks_meangs", type=int, default=16,
                   help="Max concurrent MEANGS slots  (default: %(default)s)")
    g.add_argument("--max_tasks_mt",     type=int, default=16,
                   help="Max concurrent MT slots  (default: %(default)s)")
    g.add_argument("--max_tasks_nr",     type=int, default=16,
                   help="Max concurrent NR slots  (default: %(default)s)")

    gm = sp.add_argument_group("Memory gating (optional Snakemake resource constraints)")
    gm.add_argument("--mem_total_mb",    type=int, default=None,
                    help="Total memory budget in MB.  Omit → auto-detect (85%% of system RAM on Linux)")
    gm.add_argument("--mem_meangs_mb",   type=int, default=12000, metavar="MB")
    gm.add_argument("--mem_mt_mb",       type=int, default=24000, metavar="MB")
    gm.add_argument("--mem_nr_mb",       type=int, default=24000, metavar="MB")


# ── Epilog strings ────────────────────────────────────────────────────────────
_EPILOG_MTNR = textwrap.dedent("""\
Examples
--------
  # Recommended: MT + NR in parallel, 4 samples at a time
  python3 %(prog)s -i 01.cleandata -o results \\
      --nr_seed_fasta /db/rRNA/seed.fasta \\
      --nr_genes_fasta /db/rRNA/label.fasta \\
      --max_parallel 4

  # Faster MT by providing a seed from a related species
  python3 %(prog)s -i 01.cleandata -o results \\
      --mt_seed_fasta /db/Mytilus_mt.fasta \\
      --max_reads_mt 1500000 --R_mt 5 \\
      --nr_seed_fasta /db/rRNA/seed.fasta \\
      --nr_genes_fasta /db/rRNA/label.fasta

  # Flat read layout with custom suffixes
  python3 %(prog)s -i reads_flat -o results --fq_position flat \\
      --suffix_fq "_1.clean.fq.gz,_2.clean.fq.gz;_R1.fastq.gz,_R2.fastq.gz" \\
      --nr_seed_fasta /db/rRNA/seed.fasta \\
      --nr_genes_fasta /db/rRNA/label.fasta

  # Dry-run (show plan only)
  python3 %(prog)s -i 01.cleandata -o results --dryrun \\
      --nr_seed_fasta /db/rRNA/seed.fasta --nr_genes_fasta /db/rRNA/label.fasta

Notes
-----
  • --nr_seed_fasta  and --nr_genes_fasta are REQUIRED for the NR step.
    Typical preparation:
      1. Download 18S + 28S sequences from NCBI for your taxon.
      2. Concatenate into seed.fasta (--nr_seed_fasta).
      3. Annotate/label them → label.fasta (--nr_genes_fasta).
  • MT and NR jobs run independently in parallel (no interdependency).
""")

_EPILOG_MT = textwrap.dedent("""\
Examples
--------
  # Basic MT assembly (uses built-in GetOrganelle animal_mt DB)
  python3 %(prog)s -i 01.cleandata -o mt_out

  # Faster: cap reads and provide a seed FASTA
  python3 %(prog)s -i 01.cleandata -o mt_out \\
      --mt_seed_fasta /db/close_relative_mt.fasta \\
      --max_reads_mt 1500000 --R_mt 5 --t_mt 16

  # Non-animal: plant plastid
  python3 %(prog)s -i 01.cleandata -o cp_out --F_mt embplant_pt

GetOrganelle -s (seed) speed tips
----------------------------------
  Providing a seed sequence from a closely related species via --mt_seed_fasta
  substantially speeds up read filtering because GetOrganelle can start extending
  from matching reads immediately.

  The built-in DB (animal_mt) is used when no seed is given and is sufficient
  for most metazoans, but a species-specific seed reduces run time significantly.

  Recommended:
    • --mt_seed_fasta  + --max_reads_mt 1500000  for typical 150bp WGS (>20x MT)
    • --w_mt 105       for larger word size (faster, slightly less sensitive)
    • --R_mt 5         reduces extension rounds (10 is conservative)
""")

_EPILOG_NR = textwrap.dedent("""\
Examples
--------
  # NR/rRNA assembly (18S, 28S, ITS)
  python3 %(prog)s -i 01.cleandata -o nr_out \\
      --nr_seed_fasta /db/rRNA/seed.fasta \\
      --nr_genes_fasta /db/rRNA/label.fasta

  # With read cap and extension limit (recommended for 18S/28S ~1800 bp)
  python3 %(prog)s -i 01.cleandata -o nr_out \\
      --nr_seed_fasta /db/rRNA/seed.fasta \\
      --nr_genes_fasta /db/rRNA/label.fasta \\
      --max_reads_nr 2000000 \\
      --max_extending_len_nr 250 \\
      --P_nr 0

GetOrganelle -F anonym requirements
-------------------------------------
  -F anonym REQUIRES both -s (--nr_seed_fasta) and --genes (--nr_genes_fasta).
  This follows the GetOrganelle FAQ:
    https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#custom-loci

  Prepare databases:
    1.  Collect reference rRNA sequences (18S / 28S / ITS) from NCBI.
    2.  Save as seed.fasta   → --nr_seed_fasta
    3.  Save as label.fasta  → --nr_genes_fasta
        (same sequences are acceptable for -s and --genes)
""")

_EPILOG_ALL = textwrap.dedent("""\
Examples
--------
  # Full pipeline: MEANGS -> MT -> NR
  python3 %(prog)s -i 01.cleandata -o out_all \\
      --nr_seed_fasta /db/rRNA/seed.fasta \\
      --nr_genes_fasta /db/rRNA/label.fasta \\
      --max_parallel 2

  # With MEANGS seed used for MT (default ON)
  python3 %(prog)s -i 01.cleandata -o out_all \\
      --use_meangs_seed \\
      --nr_seed_fasta /db/rRNA/seed.fasta \\
      --nr_genes_fasta /db/rRNA/label.fasta

Pipeline order (per sample)
----------------------------
  MEANGS ──► MT (may use MEANGS scaffold_seeds as -s)
         └──► NR (requires --nr_seed_fasta + --nr_genes_fasta)
  Note: MT and NR for each sample run sequentially to manage resources.
  Different samples run in parallel (governed by --max_parallel).
""")

_EPILOG_MGS2GET = textwrap.dedent("""\
Examples
--------
  # MEANGS -> GetOrganelle MT
  python3 %(prog)s -i 01.cleandata -o out

  # With MEANGS seed forwarded to GetOrganelle (default ON)
  python3 %(prog)s -i 01.cleandata -o out --use_meangs_seed

  # Faster MT via read cap
  python3 %(prog)s -i 01.cleandata -o out --max_reads_mt 1500000 --R_mt 5
""")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="meangs2get_snakemake_v4.py",
        description=textwrap.dedent("""\
            Snakemake batch pipeline for organelle / rRNA assembly from Illumina reads.

            Subcommands (modes):
              mt-nr        GetOrganelle MT + NR/rRNA in parallel          [RECOMMENDED]
              mgs2get      MEANGS -> GetOrganelle MT
              getorganelle GetOrganelle MT only
              nr-only      GetOrganelle NR/rRNA only  (-F anonym)
              all          MEANGS -> MT -> NR
            """),
        formatter_class=_HelpFormatter,
        epilog=textwrap.dedent("""\
            Run  python3 %(prog)s <subcommand> -h  for full per-mode help and examples.
            """),
    )

    # Hidden internal dispatch arguments (used by Snakemake shell rules)
    p.add_argument("--smk_internal",
                   choices=["meangs_one", "mt_one", "nr_one",
                             "finalize_mt", "finalize_nr", "finalize_mt_nr"],
                   help=argparse.SUPPRESS)
    p.add_argument("--config",  help=argparse.SUPPRESS)
    p.add_argument("--sample",  help=argparse.SUPPRESS)

    sub = p.add_subparsers(dest="cmd", metavar="<subcommand>")

    # ── mt-nr ─────────────────────────────────────────────────────────────
    mn = sub.add_parser(
        "mt-nr",
        help="GetOrganelle MT + NR/rRNA from reads  [RECOMMENDED]",
        description="Run GetOrganelle for both mitochondrial (MT) and nuclear rRNA (NR) "
                    "assemblies in parallel from paired Illumina reads.",
        epilog=_EPILOG_MTNR,
        formatter_class=_HelpFormatter,
    )
    _add_common_io(mn); _add_reads_layout(mn); _add_exec_control(mn)
    _add_mt_params(mn); _add_nr_params(mn); _add_outputs(mn); _add_smk_sched(mn)

    # ── mgs2get ───────────────────────────────────────────────────────────
    m = sub.add_parser(
        "mgs2get",
        help="MEANGS -> GetOrganelle MT",
        description="Run MEANGS first (isolated per-sample), then use the MEANGS "
                    "scaffold_seeds as a seed for GetOrganelle MT assembly.",
        epilog=_EPILOG_MGS2GET,
        formatter_class=_HelpFormatter,
    )
    _add_common_io(m); _add_reads_layout(m); _add_exec_control(m)
    _add_meangs_params(m); _add_mt_params(m); _add_outputs(m); _add_smk_sched(m)

    # ── getorganelle (MT only) ─────────────────────────────────────────────
    gt = sub.add_parser(
        "getorganelle",
        help="GetOrganelle MT only (no MEANGS)",
        description="Run GetOrganelle for mitochondrial genome assembly only.",
        epilog=_EPILOG_MT,
        formatter_class=_HelpFormatter,
    )
    _add_common_io(gt); _add_reads_layout(gt); _add_exec_control(gt)
    _add_mt_params(gt); _add_outputs(gt); _add_smk_sched(gt)

    # ── nr-only ───────────────────────────────────────────────────────────
    nr = sub.add_parser(
        "nr-only",
        help="GetOrganelle NR/rRNA only  (-F anonym + -s + --genes)",
        description="Assemble nuclear rRNA loci (18S, 28S, ITS, etc.) using "
                    "GetOrganelle with custom seed and gene databases.",
        epilog=_EPILOG_NR,
        formatter_class=_HelpFormatter,
    )
    _add_common_io(nr); _add_reads_layout(nr); _add_exec_control(nr)
    _add_nr_params(nr); _add_outputs(nr); _add_smk_sched(nr)

    # ── all ───────────────────────────────────────────────────────────────
    a = sub.add_parser(
        "all",
        help="MEANGS -> GetOrganelle MT -> GetOrganelle NR (full pipeline)",
        description="Full pipeline: MEANGS, then MT assembly (optionally seeded from MEANGS), "
                    "then NR/rRNA assembly.",
        epilog=_EPILOG_ALL,
        formatter_class=_HelpFormatter,
    )
    _add_common_io(a); _add_reads_layout(a); _add_exec_control(a)
    _add_meangs_params(a); _add_mt_params(a); _add_nr_params(a)
    _add_outputs(a); _add_smk_sched(a)

    return p


# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
def _preflight_tools(args, need_meangs: bool) -> None:
    if need_meangs:
        build_tool_cmd(
            tool_label="MEANGS",
            exe_name=MEANGS_EXE_NAME,
            base_args=["-h"],
            conda_env=getattr(args, "conda_env_meangs", None),
            tool_path=getattr(args, "meangs_path", None),
            conda_exe=getattr(args, "conda_exe", "conda"),
        )
    build_tool_cmd(
        tool_label="GetOrganelle",
        exe_name=GETORG_READS_EXE_NAME,
        base_args=["-h"],
        conda_env=getattr(args, "conda_env_getorganelle", None),
        tool_path=getattr(args, "getorganelle_path", None),
        conda_exe=getattr(args, "conda_exe", "conda"),
    )


# ---------------------------------------------------------------------------
# Sample collection
# ---------------------------------------------------------------------------
def _collect_samples(
    args,
) -> Tuple[Dict[str, dict], List[str], List[Tuple[str, Optional[Path], Optional[Path]]]]:
    input_dir = Path(args.input_dir).resolve()

    # Parse suffix pairs (new multi-suffix support)
    try:
        suffix_pairs = parse_suffix_pairs(args.suffix_fq)
    except ValueError as exc:
        logging.error(f"[SUFFIX] {exc}")
        raise

    logging.info(
        f"[SUFFIX] using {len(suffix_pairs)} pair(s): "
        + " | ".join(f"{r1} + {r2}" for r1, r2 in suffix_pairs)
    )

    subdir_list = read_list_maybe_txt(args.subdir_list, input_dir)
    tasks = collect_tasks(input_dir, args.fq_position, subdir_list, suffix_pairs)
    logging.info(f"[TASKS] collected={len(tasks)} fq_position={args.fq_position}")
    if not tasks:
        return {}, [], tasks

    samples: Dict[str, dict] = {}
    missing_reads: List[str] = []
    for tid, fq1, fq2 in tasks:
        sid = safe_name(tid)
        if fq1 is None or fq2 is None:
            logging.warning(f"[MISSING_READS] {sid}: no paired reads found")
            missing_reads.append(sid)
        else:
            samples[sid] = {
                "r1": str(Path(fq1).resolve()),
                "r2": str(Path(fq2).resolve()),
            }
    return samples, missing_reads, tasks


# ---------------------------------------------------------------------------
# Config assembly
# ---------------------------------------------------------------------------
def _make_common_params(args, base_dir_for_paths: Path) -> dict:
    nr_seed  = resolve_path_maybe_relative(getattr(args, "nr_seed_fasta",  None), base_dir_for_paths)
    nr_genes = resolve_path_maybe_relative(getattr(args, "nr_genes_fasta", None), base_dir_for_paths)
    mt_seed  = resolve_path_maybe_relative(getattr(args, "mt_seed_fasta",  None), base_dir_for_paths)

    return {
        "conda_exe":               args.conda_exe,
        "conda_env_meangs":        getattr(args, "conda_env_meangs", None),
        "conda_env_getorganelle":  args.conda_env_getorganelle,
        "meangs_path":             getattr(args, "meangs_path", None),
        "getorganelle_path":       args.getorganelle_path,

        "clean_get":    bool(args.clean_get),
        "clean_meangs": bool(getattr(args, "clean_meangs", True)),

        # MT parameters
        "F_mt":          getattr(args, "F_mt",  "animal_mt"),
        "R_mt":          int(getattr(args, "R_mt",  10)),
        "t_mt":          int(getattr(args, "t_mt",   8)),
        "k_mt":          getattr(args, "k_mt",  "21,45,65,85,105"),
        "w_mt":          getattr(args, "w_mt",  None),
        "max_reads_mt":  getattr(args, "max_reads_mt", None),
        "all_data_mt":   bool(getattr(args, "all_data_mt", False)),
        "mt_seed_fasta": mt_seed,
        "get_extra_mt":  list(getattr(args, "get_extra_mt", []) or []),

        "use_meangs_seed":   bool(getattr(args, "use_meangs_seed", False)),
        "seed_glob_meangs":  getattr(args, "seed_glob_meangs", None),

        # NR parameters
        "F_nr":                  getattr(args, "F_nr",  "anonym"),
        "nr_seed_fasta":         nr_seed,
        "nr_genes_fasta":        nr_genes,
        "R_nr":                  int(getattr(args, "R_nr",  10)),
        "t_nr":                  int(getattr(args, "t_nr",   8)),
        "k_nr":                  getattr(args, "k_nr",  "21,45,65,85,105"),
        "P_nr":                  getattr(args, "P_nr",  "0"),
        "max_extending_len_nr":  getattr(args, "max_extending_len_nr", ""),
        "max_reads_nr":          getattr(args, "max_reads_nr", None),
        "all_data_nr":           bool(getattr(args, "all_data_nr", False)),
        "get_extra_nr":          list(getattr(args, "get_extra_nr", []) or []),

        # Memory
        "mem_meangs_mb": int(getattr(args, "mem_meangs_mb", 12000)),
        "mem_mt_mb":     int(getattr(args, "mem_mt_mb",     24000)),
        "mem_nr_mb":     int(getattr(args, "mem_nr_mb",     24000)),
    }


def _make_cfg(
    out_root: Path,
    args,
    samples: Dict[str, dict],
    missing_reads: List[str],
    mode: str,
    use_meangs_for_mt: bool,
) -> dict:
    return {
        "driver_py":       str(Path(__file__).resolve()),
        "python_exe":      sys.executable,
        "out_root":        str(out_root),
        "meangs_dir":      args.meangs_dir,
        "mt_dir":          args.mt_dir,
        "nr_dir":          args.nr_dir,
        "samples":         samples,
        "missing_reads":   sorted(set(missing_reads)),
        "use_meangs_for_mt": bool(use_meangs_for_mt),

        "params": _make_common_params(args, base_dir_for_paths=Path.cwd()),

        "summary_mt":  str((out_root / args.summary_mt_name).resolve()),
        "summary_nr":  str((out_root / args.summary_nr_name).resolve()),
        "summary_all": str((out_root / args.summary_all_name).resolve()),
        "report_tsv":  str((out_root / args.report_name).resolve()),

        "failed_mt":   str((out_root / args.failed_mt_name).resolve()),
        "failed_nr":   str((out_root / args.failed_nr_name).resolve()),
        "failed_all":  str((out_root / args.failed_all_name).resolve()),

        "mode": mode,
    }


def _compute_cores(args) -> int:
    auto = os.cpu_count() or 8
    cores = int(args.smk_cores) if int(args.smk_cores) > 0 else auto
    max_rule = max(
        int(getattr(args, "threads_meangs", 1)),
        int(getattr(args, "t_mt", 1)),
        int(getattr(args, "t_nr", 1)),
    )
    return max(cores, max_rule)


def _auto_mem_total(args) -> Optional[int]:
    m = getattr(args, "mem_total_mb", None)
    if m is not None:
        return int(m)
    auto_mb = detect_total_mem_mb_linux()
    if auto_mb:
        return int(auto_mb * 0.85)
    return None


# ---------------------------------------------------------------------------
# Main workflow runner
# ---------------------------------------------------------------------------
def _run_mode(args, mode: str, need_meangs: bool, use_meangs_for_mt: bool) -> int:
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_logger(out_root / "pipeline.log")

    # Tool pre-flight
    try:
        _preflight_tools(args, need_meangs=need_meangs)
    except FileNotFoundError as e:
        logging.error(str(e))
        return 2

    # Collect sample reads
    try:
        samples, missing_reads, tasks = _collect_samples(args)
    except ValueError:
        return 2

    if not tasks:
        logging.error("[TASKS] no tasks found in input directory.")
        return 1
    if not samples:
        logging.error("[TASKS] no valid paired reads found; all samples missing reads.")
        write_list(sorted(set(missing_reads)), out_root / args.failed_all_name)
        return 1

    # Validate NR databases if needed
    if mode in ("nr-only", "mt-nr", "all"):
        seed  = resolve_path_maybe_relative(getattr(args, "nr_seed_fasta",  None), Path.cwd())
        genes = resolve_path_maybe_relative(getattr(args, "nr_genes_fasta", None), Path.cwd())
        if not seed or not Path(seed).is_file():
            logging.error(
                "[NR] --nr_seed_fasta is REQUIRED for NR modes and must point to an existing file.\n"
                "     Typical usage: --nr_seed_fasta /db/rRNA/seed.fasta"
            )
            return 2
        if not genes or not Path(genes).is_file():
            logging.error(
                "[NR] --nr_genes_fasta is REQUIRED for NR modes and must point to an existing file.\n"
                "     Typical usage: --nr_genes_fasta /db/rRNA/label.fasta"
            )
            return 2

    cfg = _make_cfg(out_root, args, samples, missing_reads,
                    mode=mode, use_meangs_for_mt=use_meangs_for_mt)

    # Add MEANGS-specific params to config when needed
    if need_meangs:
        cfg["params"].update({
            "meangs_outname":       args.meangs_outname,
            "threads_meangs":       int(args.threads_meangs),
            "insert_meangs":        int(args.insert_meangs),
            "species_class_meangs": args.species_class_meangs,
            "nsample_meangs":       int(args.nsample_meangs),
            "seqscaf_meangs":       args.seqscaf_meangs,
            "keepMinLen_meangs":    args.keepMinLen_meangs,
            "clip_meangs":          bool(args.clip_meangs),
            "keepIntMed_meangs":    bool(args.keepIntMed_meangs),
            "skipassem_meangs":     bool(args.skipassem_meangs),
            "skipqc_meangs":        bool(args.skipqc_meangs),
            "skiphmm_meangs":       bool(args.skiphmm_meangs),
            "skipextend_meangs":    bool(args.skipextend_meangs),
            "deepin_meangs":        bool(args.deepin_meangs),
            "meangs_extra":         list(args.meangs_extra or []),
        })

    # Map subcommand mode to Snakefile mode string
    snake_mode = {
        "mgs2get":      "mgs2get",
        "getorganelle": "mt-only",
        "nr-only":      "nr-only",
        "mt-nr":        "mt-nr",
        "all":          "all",
    }[mode]

    workflow_dir = out_root / ".smk"
    snakefile, _ = write_smk_files(workflow_dir, cfg, mode=snake_mode)

    logging.info(f"[SMK] workflow_dir={workflow_dir}")
    logging.info(f"[SMK] mode={mode}  samples={len(samples)}  missing_reads={len(missing_reads)}")

    cores        = _compute_cores(args)
    mem_total_mb = _auto_mem_total(args)
    if mem_total_mb:
        logging.info(f"[SCHED] mem_total_mb={mem_total_mb}")

    rc = run_snakemake_cli(
        workflow_dir=workflow_dir,
        snakefile=snakefile,
        cores=cores,
        jobs=int(args.max_parallel),
        latency_wait=int(args.latency_wait),
        dryrun=bool(args.dryrun),
        mem_total_mb=mem_total_mb,
        max_tasks_meangs=int(args.max_tasks_meangs),
        max_tasks_mt=int(args.max_tasks_mt),
        max_tasks_nr=int(args.max_tasks_nr),
    )

    if rc == 0 and not args.dryrun:
        failed_path = out_root / args.failed_all_name
        n_failed = 0
        if failed_path.exists():
            n_failed = len([
                x for x in failed_path.read_text(encoding="utf-8", errors="ignore").splitlines()
                if x.strip()
            ])
        logging.info(f"[DONE] Pipeline finished. failed_samples={n_failed}")
        logging.info(f"[OUT] Summary FASTA : {out_root / args.summary_all_name}")
        logging.info(f"[OUT] Report TSV    : {out_root / args.report_name}")
        if n_failed:
            logging.warning(f"[WARN] {n_failed} sample(s) failed → {failed_path}")
    return rc


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main() -> int:
    parser = build_parser()
    args   = parser.parse_args()

    # ── Internal Snakemake dispatch ────────────────────────────────────────
    if getattr(args, "smk_internal", None):
        config_path = Path(args.config).resolve() if args.config else None
        if not config_path or not config_path.exists():
            print("[FATAL] --smk_internal requires --config <config.json>", file=sys.stderr)
            return 2
        dispatch = {
            "meangs_one":     lambda: _smk_meangs_one(config_path, args.sample),
            "mt_one":         lambda: _smk_mt_one(config_path, args.sample),
            "nr_one":         lambda: _smk_nr_one(config_path, args.sample),
            "finalize_mt":    lambda: _smk_finalize_mt(config_path),
            "finalize_nr":    lambda: _smk_finalize_nr(config_path),
            "finalize_mt_nr": lambda: _smk_finalize_mt_nr(config_path),
        }
        fn = dispatch.get(args.smk_internal)
        if fn is None:
            print(f"[FATAL] unknown --smk_internal: {args.smk_internal}", file=sys.stderr)
            return 2
        return fn()

    # ── Subcommand dispatch ────────────────────────────────────────────────
    _CMD_TABLE = {
        "mgs2get":      dict(mode="mgs2get",      need_meangs=True,  use_meangs_for_mt=True),
        "getorganelle": dict(mode="getorganelle",  need_meangs=False, use_meangs_for_mt=False),
        "nr-only":      dict(mode="nr-only",       need_meangs=False, use_meangs_for_mt=False),
        "mt-nr":        dict(mode="mt-nr",         need_meangs=False, use_meangs_for_mt=False),
        "all":          dict(mode="all",           need_meangs=True,  use_meangs_for_mt=True),
    }

    if args.cmd in _CMD_TABLE:
        return _run_mode(args, **_CMD_TABLE[args.cmd])

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
