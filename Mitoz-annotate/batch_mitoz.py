#!/usr/bin/env python3
# DESC: Integrated MitoZ batch annotation pipeline
# VERSION: v4 (Added: 20260325_1556)
# -*- coding: utf-8 -*-
"""
Integrated MitoZ batch annotation pipeline — publication-grade.

v0.28 over v0.27:
  BUG-V0.28-A  _acc_variants: return type changed from List[str] to
               Tuple[str, ...] so the object stored by lru_cache is
               immutable.  The previous List was safe only by convention
               ("callers must not mutate the result"); any accidental
               append/clear/sort would silently corrupt the cache and cause
               spurious lookup failures.  All call-sites only iterate over
               the result, so this is a safe drop-in fix.
  BUG-V0.28-B  workflow_replace_organism TSV column mismatch fixed.
               The header previously declared 5 columns but exception rows
               wrote 6 (the FAIL message was appended to the numeric
               ``total`` column, breaking downstream TSV parsers).  A
               dedicated ``note`` column is now declared in both header and
               all data rows; normal rows carry an empty note field.
  BUG-V0.28-C  rewrite_genbank_record_ids: replaced ``rr = r[:]``
               slice-copy with a direct SeqRecord construction (same
               pattern as BUG-NEW-E in v0.27).  The slice triggered
               Biopython's __getitem__, which clones every feature location,
               only for the result to be discarded when .id / .name /
               .annotations were immediately overwritten.
  BUG-V0.28-D  find_any_summary_txt: sort key changed from
               ``len(str(p))`` (raw string length) to ``len(p.parts)``
               (actual directory depth).  String length is not a reliable
               proxy for depth when parent directory names differ in length.
  BUG-V0.28-F  RunSummary: added ``reorient_failed`` field.
               workflow_reorient was appending to ``annotation_failed``,
               which is semantically wrong (no MitoZ annotation step exists
               in the reorient workflow).  ``reorient_failed`` is included
               in all_problem_samples() and log().
  NOTE         BUG-V0.28-E (interactive_resolve KeyError) was assessed as
               a false positive: the existing ``exact_dict.get(chosen)``
               call already handles the case correctly.

v0.27 over v0.26:
  BUG-NEW-A  _safe_copy: promoted from a nested closure inside
             workflow_transfer_metadata to a module-level utility.  All
             shutil.copy2 call-sites that could theoretically receive the
             same src/dst path (e.g. when --out_root == input dir) are now
             protected against shutil.SameFileError uniformly.
  BUG-NEW-B  build_organism_index: removed the redundant
             ``index.setdefault(norm, entry)`` executed after the
             ``_acc_variants`` loop.  _acc_variants already registers the
             normalised form; the extra call was dead code and could mask
             later entries with the wrong value when the same norm key
             appeared under multiple accessions.
  BUG-NEW-C  _split_gb_records: filter out empty / whitespace-only
             pseudo-records produced when a GenBank file ends with trailing
             newlines after the final ``//`` terminator.  Previously those
             empty strings propagated into replace_metadata_in_text,
             incrementing stats["total_new"] and stats["no_accession"]
             by one phantom record each, causing misleading summary counts.
  BUG-NEW-D  workflow_annotate / workflow_run merge consistency: both
             workflows now use rglob("*.gbf") when building all.final.gbf
             so that no_{gene}_gene/ subdirectory records are included in
             both cases.  Previous behaviour (glob in annotate, rglob in
             run) was documented but confusing; the NOTE at the rglob call
             is retained to make the intentional include explicit.
  BUG-NEW-E  replace_organism_in_record / _sync_organism_qualifier:
             eliminated the wasteful rec[:] slice-copy whose only output
             (the feature list with original locations) was immediately
             discarded when new_rec.features = new_features was assigned.
             The SeqRecord is now built directly from its constituent parts,
             saving one full feature-copy per record.
  BUG-NEW-F  finalize_genbank_output Steps 3 / 3.5: replaced bare
             shutil.copy2 calls with the module-level _safe_copy so that
             SameFileError is silently ignored if the caller happens to
             point dest_gb at the same path as working_gb (edge case when
             out_root == annotation workdir).
  OPT-1      replace_metadata_in_text: pre-compile the LOCUS bp-count
             substitution regex once per call rather than once per record.
  OPT-2      _acc_variants: added @functools.lru_cache so repeated calls
             with the same accession string (common during index building)
             are free after the first computation.
  NOTE       Snakemake support: the --scheduler_backend argument controls
             dispatch.  Use ``--scheduler_backend snakemake`` (or ``auto``
             which picks Snakemake when available, threads otherwise).
             The embedded worker script (_mitoz_worker.py) runs MitoZ
             annotation only; finalize_genbank_output / organism-sync are
             always executed in the main process after all workers finish,
             so all v0.26 / v0.27 bug fixes apply regardless of scheduler.

v0.26 over v0.25:
  BUG-10 build_organism_index: now returns Dict[str, Tuple[str, str]]
  BUG-11 replace_organism_in_record: removed erroneous ≤8-word heuristic
  BUG-12 finalize_genbank_output: added Step 3.5 — sync /organism=
  NEW    _sync_organism_qualifier helper
  BUG-13 lookup_organism: return type Optional[Tuple[str, str]]
  BUG-14 workflow_replace_organism: updated for Tuple return values
  BUG-A  replace_organism_in_record: removed "organism" in new_quals guard
  BUG-B  _sync_organism_qualifier: same guard removed

v0.25 over v0.24:
  NEW    replace_organism subcommand

v0.24 over v0.23:
  HOLE-1 gb_to_internal_fasta: write all contigs
  HOLE-2 workflow_run Phase-1b: only blacklist when ALL contigs lack gene
  HOLE-3 split_or_sanitize_fastas: global seen-set deduplication

v0.23 over v0.22:
  BUG-8  rewrite_genbank_record_ids: force-update annotations['accessions']
  BUG-9  finalize_genbank_output Step 2: merge reoriented + un-reoriented

v0.22 over v0.21:
  BUG-6  _get_record_accession: scan VERSION > ACCESSION > LOCUS
  BUG-7  shift_location: zero-length point feature fix (start <= end)
  OPT-1  LOCUS-line bp-count regex broadened

v0.21 over v0.20:
  BUG-5a/b/c shift_location / workflow_transfer_metadata fixes

v0.20 over v0.19:
  BUG-1..4

Features
--------
1. annotate          — batch MitoZ annotate with internal short IDs
2. reorient          — reorient existing .gbf to start at a target gene
3. run               — annotate → reorient → optional re-annotation →
                       guaranteed final reorientation → metadata transfer
                       → organism/source qualifier sync (Step 3.5)
4. transfer_metadata — standalone GenBank header-transfer tool
5. replace_organism  — restore /organism= in FEATURES source from original GB

Scheduler backends
------------------
``--scheduler_backend auto``      (default) use Snakemake if available,
                                   else ThreadPoolExecutor.
``--scheduler_backend snakemake`` always use Snakemake.
``--scheduler_backend threads``   always use ThreadPoolExecutor.

The Snakemake path generates:
  _snakemake_{phase}/Snakefile
  _snakemake_{phase}/_mitoz_worker.py  (self-contained worker; annotation only)
  _snakemake_{phase}/{id}.json         (per-sample task payloads)
Finalization (ID-rewrite, reorientation, metadata overlay, organism sync)
is executed in the main process after all Snakemake jobs complete, so all
pipeline bug-fixes apply regardless of backend.

Metadata preservation guarantee (v0.27)
----------------------------------------
When GenBank files are used as input and --preserve_phase1_metadata yes
(the default for 'run'), the final output records are guaranteed to have:
  • All header fields (LOCUS bp-count, DEFINITION, ACCESSION, VERSION,
    KEYWORDS, SOURCE line, ORGANISM + taxonomy lineage, REFERENCE,
    COMMENT, DBLINK) — from the original input GenBank.
  • Re-annotated gene features (CDS / tRNA / rRNA) — from MitoZ.
  • /organism= qualifier inside FEATURES source — from the original input
    GenBank (synced in Step 3.5, fixes MitoZ placeholder strings).
  • /source= and all other source-feature qualifiers — from MitoZ.

Requirements
------------
Python ≥ 3.8, biopython, mitoz in PATH or via --conda_env_mitoz/--mitoz_path,
snakemake (optional).

Citation
--------
MitoZ:      Meng et al. (2019) NAR doi:10.1093/nar/gkz173
Snakemake:  Mölder et al. (2021) F1000Res doi:10.12688/f1000research.29032.2
Biopython:  Cock et al. (2009) Bioinformatics doi:10.1093/bioinformatics/btp163
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import difflib
import functools
import json
import logging
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
except Exception:
    SeqIO = None  # type: ignore

MITOZ_EXE_NAME = "mitoz"

# ─────────────────────────────────────────────────────────────────────────────
# Logging
# ─────────────────────────────────────────────────────────────────────────────
_LOG_LOCK = threading.Lock()


class _LockedHandler(logging.Handler):
    def __init__(self, inner: logging.Handler):
        super().__init__()
        self._inner = inner
        self.setFormatter(inner.formatter)
        self.setLevel(inner.level)

    def emit(self, record: logging.LogRecord) -> None:
        with _LOG_LOCK:
            self._inner.emit(record)


def setup_root_logger(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    for h in list(root.handlers):
        root.removeHandler(h)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh = logging.FileHandler(str(log_file), encoding="utf-8")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    root.addHandler(_LockedHandler(fh))
    root.addHandler(_LockedHandler(sh))


def get_sample_logger(name: str, log_file: Path) -> logging.Logger:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    lg = logging.getLogger(name)
    lg.setLevel(logging.INFO)
    lg.propagate = False
    for h in list(lg.handlers):
        try:
            h.close()
        except Exception:
            pass
        lg.removeHandler(h)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh = logging.FileHandler(str(log_file), encoding="utf-8")
    fh.setFormatter(fmt)
    lg.addHandler(fh)
    return lg


def close_sample_logger(name: str) -> None:
    lg = logging.getLogger(name)
    for h in list(lg.handlers):
        try:
            h.close()
        except Exception:
            pass
        lg.removeHandler(h)


# ─────────────────────────────────────────────────────────────────────────────
# Generic helpers
# ─────────────────────────────────────────────────────────────────────────────
def tail_text(s: str, n_lines: int = 120) -> str:
    return "\n".join((s or "").splitlines()[-n_lines:])


def safe_name(s: str) -> str:
    s = str(s).strip().replace(" ", "_")
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", s)
    return s or "NA"


def header_first_token(desc_or_id: str) -> str:
    s = (desc_or_id or "").strip()
    if s.startswith(">"):
        s = s[1:]
    return s.split()[0] if s else "NA"


def which_or_none(exe: str) -> Optional[str]:
    try:
        return shutil.which(exe)
    except Exception:
        return None


def resolve_exe_from_path(exe_name: str, tool_path: Optional[str]) -> Optional[Path]:
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
    if conda_env:
        conda_bin = which_or_none(conda_exe)
        if not conda_bin:
            raise FileNotFoundError(f"[{tool_label}] '{conda_exe}' not found in PATH.")
        return [conda_bin, "run", "-n", conda_env, "--no-capture-output", exe_name] + base_args
    exe_path = resolve_exe_from_path(exe_name, tool_path)
    if exe_path is not None:
        if not exe_path.exists():
            raise FileNotFoundError(f"[{tool_label}] not found at: {exe_path}")
        if exe_path.suffix.lower() == ".py":
            return [sys.executable, str(exe_path)] + base_args
        return [str(exe_path)] + base_args
    found = which_or_none(exe_name)
    if not found:
        raise FileNotFoundError(
            f"[{tool_label}] cannot find '{exe_name}' in PATH.\n"
            f"  Use --conda_env_mitoz ENV  or  --mitoz_path DIR_OR_FILE"
        )
    return [found] + base_args


_FASTA_EXTS = {".fa", ".fasta", ".fna", ".fas"}
_GB_EXTS    = {".gb", ".gbk", ".genbank", ".gbf"}


def _is_fasta_file(p: Path) -> bool:
    return p.suffix.lower() in _FASTA_EXTS


def _is_gb_file(p: Path) -> bool:
    return p.suffix.lower() in _GB_EXTS


def _sniff_input_type(in_path: Path) -> str:
    if in_path.is_file():
        if _is_fasta_file(in_path):
            return "fasta"
        if _is_gb_file(in_path):
            return "genbank"
        try:
            with in_path.open(encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    line = line.strip()
                    if line:
                        if line.startswith(">"):
                            return "fasta"
                        if line.startswith("LOCUS"):
                            return "genbank"
                        break
        except Exception:
            pass
        return "fasta"
    if in_path.is_dir():
        has_fasta = any(_is_fasta_file(p) for p in in_path.rglob("*") if p.is_file())
        has_gb    = any(_is_gb_file(p)    for p in in_path.rglob("*") if p.is_file())
        if has_gb and not has_fasta:
            return "genbank"
        if has_fasta and not has_gb:
            return "fasta"
        if has_gb and has_fasta:
            return "mixed"
        return "fasta"
    raise FileNotFoundError(f"Input not found: {in_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Output helpers — only write non-empty content
# ─────────────────────────────────────────────────────────────────────────────
def write_nonempty_list(out_path: Path, items: List[str], comment: str = "") -> None:
    unique = sorted(set(items))
    if not unique:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    if comment:
        lines.append(comment + "\n")
    for item in unique:
        lines.append(item + "\n")
    out_path.write_text("".join(lines), encoding="utf-8")


def write_nonempty_table(out_path: Path, header: str, rows: List[str]) -> None:
    if not rows:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        fh.write(header + "\n")
        for r in rows:
            fh.write(r + "\n")


def int_mem_mb(x, minimum: int = 128) -> int:
    try:
        v = int(math.ceil(float(x)))
    except Exception:
        v = minimum
    return max(minimum, v)


def positive_int_arg(x: str) -> int:
    try:
        v = int(x)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(f"expected positive integer, got {x!r}") from exc
    if v < 1:
        raise argparse.ArgumentTypeError(f"expected positive integer, got {x!r}")
    return v


def id_maxlen_arg(x: str) -> int:
    v = positive_int_arg(x)
    if v < len("S00001"):
        raise argparse.ArgumentTypeError(
            f"--id_maxlen must be >= {len('S00001')} to keep internal IDs unique")
    return v


# ─────────────────────────────────────────────────────────────────────────────
# Module-level _safe_copy  (BUG-NEW-A v0.27)
# ─────────────────────────────────────────────────────────────────────────────
def _safe_copy(src: Path, dst: Path) -> None:
    """Copy *src* → *dst*, silently ignoring shutil.SameFileError.

    BUG-NEW-A (v0.27): promoted from a nested closure in
    workflow_transfer_metadata to a module-level utility so that all
    shutil.copy2 call-sites gain SameFileError protection uniformly.
    """
    try:
        shutil.copy2(str(src), str(dst))
    except shutil.SameFileError:
        pass


# ─────────────────────────────────────────────────────────────────────────────
# GenBank record accession helpers
# ─────────────────────────────────────────────────────────────────────────────
def gb_record_accession(rec: "SeqRecord") -> str:
    ver  = rec.annotations.get("sequence_version", "")
    accs = rec.annotations.get("accessions", [])
    acc  = accs[0] if accs else ""
    if not acc:
        acc = rec.name if rec.name and rec.name not in (".", "<unknown>") else ""
    if not acc:
        acc = rec.id  if rec.id  and rec.id  not in (".", "<unknown>") else ""
    if ver and acc and not acc.endswith(f".{ver}"):
        versioned = f"{acc}.{ver}"
    else:
        versioned = acc
    return safe_name(versioned or rec.id or "unknown")


@dataclass
class GbRecordMeta:
    acc_versioned: str
    acc_base:      str
    version_int:   int
    seq_hash:      str
    seq_len:       int
    topology:      str
    source_file:   Path
    record:        object


def _gb_version_int(rec: "SeqRecord") -> int:
    ver = rec.annotations.get("sequence_version", "")
    try:
        return int(ver)
    except (TypeError, ValueError):
        pass
    m = re.search(r"\.(\d+)$", rec.id or "")
    return int(m.group(1)) if m else 0


def dedup_gb_records(
    gb_files: List[Path],
    dedup_report_path: Optional[Path] = None,
) -> Tuple[List[GbRecordMeta], Dict[str, Path], Dict]:
    import hashlib

    stats: Dict = {
        "total_in": 0, "kept": 0,
        "dup_exact": 0, "dup_version": 0, "dup_seq": 0,
        "parse_errors": 0, "no_sequence": 0,
        "dup_exact_list": [], "dup_version_list": [], "dup_seq_list": [],
        "no_sequence_list": [],
    }

    all_metas:               List[GbRecordMeta] = []
    metadata_only_acc_to_gb: Dict[str, Path]    = {}

    for gb_path in gb_files:
        try:
            records = list(SeqIO.parse(str(gb_path), "genbank"))
        except Exception as exc:
            logging.warning(f"[DEDUP] Cannot parse {gb_path}: {exc}")
            stats["parse_errors"] += 1
            continue

        for rec in records:
            stats["total_in"] += 1
            acc_ver  = gb_record_accession(rec)
            acc_base = _normalize_acc(acc_ver) or acc_ver
            ver_int  = _gb_version_int(rec)

            try:
                seq_str = str(rec.seq).upper().replace(" ", "")
                has_seq = bool(seq_str)
            except Exception:
                seq_str, has_seq = "", False

            if not has_seq:
                stats["no_sequence"] += 1
                stats["no_sequence_list"].append(f"{acc_ver} ({gb_path.name})")
                metadata_only_acc_to_gb.setdefault(safe_name(acc_ver), gb_path.resolve())
                continue

            seq_hash = hashlib.md5(seq_str.encode()).hexdigest()
            topo     = rec.annotations.get("topology", "circular").lower()
            if topo not in ("circular", "linear"):
                topo = "circular"

            all_metas.append(GbRecordMeta(
                acc_versioned=acc_ver, acc_base=acc_base,
                version_int=ver_int, seq_hash=seq_hash,
                seq_len=len(seq_str), topology=topo,
                source_file=gb_path, record=rec,
            ))

    by_acc: Dict[str, GbRecordMeta] = {}
    for meta in all_metas:
        key = meta.acc_base
        if key not in by_acc:
            by_acc[key] = meta
        else:
            prev = by_acc[key]
            if meta.version_int > prev.version_int:
                kind = "dup_exact" if prev.acc_versioned == meta.acc_versioned else "dup_version"
                stats[kind] += 1
                stats[f"{kind}_list"].append(f"{prev.acc_versioned} superseded by {meta.acc_versioned}")
                by_acc[key] = meta
            else:
                stats["dup_exact"] += 1
                stats["dup_exact_list"].append(
                    f"{meta.acc_versioned} (dup of {prev.acc_versioned})")

    seq_hash_to_acc: Dict[str, str] = {}
    kept_records: List[GbRecordMeta] = []
    for meta in by_acc.values():
        if meta.seq_hash in seq_hash_to_acc:
            stats["dup_seq"] += 1
            stats["dup_seq_list"].append(
                f"{meta.acc_versioned} (seq identical to {seq_hash_to_acc[meta.seq_hash]}, "
                f"len={meta.seq_len}, dropped)")
            logging.warning(
                f"[DEDUP_SEQ] '{meta.acc_versioned}' has identical sequence to "
                f"'{seq_hash_to_acc[meta.seq_hash]}' (len={meta.seq_len}). Keeping only first.")
        else:
            seq_hash_to_acc[meta.seq_hash] = meta.acc_versioned
            kept_records.append(meta)

    stats["kept"] = len(kept_records)

    if dedup_report_path is not None and (kept_records or stats["total_in"] > 0):
        dedup_report_path.parent.mkdir(parents=True, exist_ok=True)
        lines = [
            "# GB deduplication report\n",
            f"# total_in={stats['total_in']}  with_seq={stats['kept']}  "
            f"no_seq(metadata_only)={stats['no_sequence']}  "
            f"dropped_dup={stats['dup_exact'] + stats['dup_version'] + stats['dup_seq']}  "
            f"parse_errors={stats['parse_errors']}\n",
            "accession\tstatus\tsource_file\tseq_len\ttopology\n",
        ]
        for meta in kept_records:
            lines.append(
                f"{meta.acc_versioned}\tkept\t{meta.source_file.name}\t"
                f"{meta.seq_len}\t{meta.topology}\n")
        for msg in stats["dup_seq_list"]:
            lines.append(f"\t\tdup_seq: {msg}\n")
        dedup_report_path.write_text("".join(lines), encoding="utf-8")

    return kept_records, metadata_only_acc_to_gb, stats


def convert_gb_to_fasta_for_annotation(
    gb_files: List[Path],
    out_dir: Path,
    dedup_report_path: Optional[Path] = None,
) -> Tuple[List[Path], Dict[str, Path]]:
    if SeqIO is None:
        raise RuntimeError("Biopython not installed: pip install biopython")
    out_dir.mkdir(parents=True, exist_ok=True)

    kept_records, metadata_only_acc_to_gb, stats = dedup_gb_records(
        gb_files, dedup_report_path=dedup_report_path)

    n_total  = stats["total_in"]
    n_kept   = stats["kept"]
    n_no_seq = stats["no_sequence"]
    n_dup    = n_total - n_kept - n_no_seq

    logging.info(
        f"[GB_DEDUP] total={n_total}  with_seq={n_kept}  "
        f"no_seq(metadata_only)={n_no_seq}  dropped_dup={n_dup}  "
        f"(dup_exact={stats['dup_exact']}, dup_version={stats['dup_version']}, "
        f"dup_seq={stats['dup_seq']}, parse_errors={stats['parse_errors']})"
    )
    if n_no_seq > 0:
        logging.info(
            f"[GB_NO_SEQ] {n_no_seq} record(s) have no sequence (ORIGIN absent) "
            f"— skipped for annotation, metadata source still tracked.")
    for msg in stats["dup_seq_list"]:
        logging.warning(f"  [DUP_SEQ]     {msg}")

    fasta_paths:    List[Path]      = []
    acc_to_orig_gb: Dict[str, Path] = {}

    for meta in kept_records:
        fa_name = safe_name(meta.acc_versioned)
        out_fa  = out_dir / f"{fa_name}.fasta"
        out_rec = SeqRecord(
            meta.record.seq, id=fa_name, name=fa_name,
            description=f"topology={meta.topology}",
        )
        with out_fa.open("w", encoding="utf-8") as fh:
            SeqIO.write([out_rec], fh, "fasta")
        fasta_paths.append(out_fa.resolve())
        acc_to_orig_gb[fa_name] = meta.source_file.resolve()

    for acc_key, src_path in metadata_only_acc_to_gb.items():
        acc_to_orig_gb.setdefault(acc_key, src_path)

    return fasta_paths, acc_to_orig_gb


def iter_inputs(
    in_path: Path,
    gb_to_fasta_dir: Optional[Path] = None,
) -> Tuple[List[Path], str, Dict[str, Path]]:
    itype = _sniff_input_type(in_path)
    logging.info(f"[INPUT] Detected input type: {itype}  ({in_path})")

    if itype == "fasta":
        if in_path.is_file():
            return [in_path.resolve()], "fasta", {}
        files = sorted(p.resolve() for p in in_path.rglob("*")
                       if p.is_file() and _is_fasta_file(p))
        if not files:
            raise FileNotFoundError(f"No FASTA files in: {in_path}")
        return files, "fasta", {}

    if itype in ("genbank", "mixed"):
        if in_path.is_file():
            gb_files = [in_path.resolve()]
        else:
            gb_files = sorted(p.resolve() for p in in_path.rglob("*")
                              if p.is_file() and _is_gb_file(p))
        if not gb_files:
            raise FileNotFoundError(f"No GenBank files in: {in_path}")
        if itype == "mixed":
            logging.warning(
                "[INPUT] Mixed FASTA+GB directory — only GB files will be used.")
        if gb_to_fasta_dir is None:
            raise RuntimeError("gb_to_fasta_dir required when input is GenBank")
        dedup_report = gb_to_fasta_dir / "dedup_report.tsv"
        fasta_paths, acc_to_orig_gb = convert_gb_to_fasta_for_annotation(
            gb_files, gb_to_fasta_dir, dedup_report_path=dedup_report)
        if not fasta_paths:
            raise RuntimeError(f"No records extracted from GenBank input: {in_path}")
        logging.info(
            f"[INPUT] Converted {len(gb_files)} GB file(s) -> "
            f"{len(fasta_paths)} unique FASTA record(s)")
        return fasta_paths, "converted_from_gb", acc_to_orig_gb

    raise FileNotFoundError(f"Input not found or unrecognised: {in_path}")


def infer_topology_from_description(desc: str) -> str:
    m = re.search(r"topology\s*=\s*(circular|linear)", (desc or "").lower())
    return m.group(1) if m else "circular"


def _is_topology_tag(s: str) -> bool:
    return bool(re.match(r"^topology\s*=\s*(circular|linear)$", (s or "").strip(), re.I))


def make_internal_id(idx: int, max_len: int) -> str:
    base = f"S{idx:05d}"
    return base[:max_len] if len(base) > max_len else base


# ─────────────────────────────────────────────────────────────────────────────
# SampleInfo / RunSummary
# ─────────────────────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class SampleInfo:
    sample_name:    str
    orig_token:     str
    internal_id:    str
    topology:       str
    source_fasta:   Path
    internal_fasta: Path
    orig_gb_path:   Optional[Path] = None


@dataclass
class RunSummary:
    annotation_failed:        List[str] = field(default_factory=list)
    final_write_failed:       List[str] = field(default_factory=list)
    no_target_gene:           List[str] = field(default_factory=list)
    metadata_transferred:     List[str] = field(default_factory=list)
    metadata_compat_match:    List[str] = field(default_factory=list)
    metadata_key_mismatch:    List[str] = field(default_factory=list)
    metadata_skipped:         List[str] = field(default_factory=list)
    metadata_transfer_failed: List[str] = field(default_factory=list)
    # BUG-V0.28-F: dedicated field for reorientation failures so that
    # workflow_reorient does not misuse annotation_failed (which implies a
    # MitoZ annotation step, absent in the reorient workflow).
    reorient_failed:          List[str] = field(default_factory=list)

    def all_problem_samples(self) -> List[str]:
        return sorted(set(
            self.annotation_failed + self.final_write_failed
            + self.no_target_gene  + self.metadata_transfer_failed
            # BUG-V0.29-I: compat/key-mismatch metadata transfers are audit
            # problems too; the old list omitted them from all_problem_samples.
            + self.metadata_key_mismatch
            + self.reorient_failed
        ))

    def log(self) -> None:
        logging.info(
            "[SUMMARY] annotation_failed=%d  final_write_failed=%d  no_target_gene=%d"
            "  reorient_failed=%d",
            len(self.annotation_failed), len(self.final_write_failed), len(self.no_target_gene),
            len(self.reorient_failed))
        logging.info(
            "[SUMMARY] metadata_transferred=%d  compat_match=%d  key_mismatch=%d  "
            "skipped=%d  transfer_failed=%d",
            len(self.metadata_transferred), len(self.metadata_compat_match),
            len(self.metadata_key_mismatch), len(self.metadata_skipped),
            len(self.metadata_transfer_failed))


# ─────────────────────────────────────────────────────────────────────────────
# FASTA helpers
# ─────────────────────────────────────────────────────────────────────────────
def write_one_record_fasta(record: "SeqRecord", out_fa: Path, seq_id: str, topology: str) -> None:
    rec = record[:]
    rec.id = rec.name = seq_id
    rec.description = f"topology={topology}"
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with out_fa.open("w", encoding="utf-8") as fh:
        SeqIO.write([rec], fh, "fasta")


def split_or_sanitize_fastas(
    in_files: List[Path],
    out_dir: Path,
    if_split: bool,
    id_maxlen: int,
    acc_to_orig_gb: Optional[Dict[str, Path]] = None,
) -> Tuple[List[SampleInfo], Path]:
    if SeqIO is None:
        raise RuntimeError("Biopython not installed: pip install biopython")
    out_dir.mkdir(parents=True, exist_ok=True)
    if acc_to_orig_gb is None:
        acc_to_orig_gb = {}
    sample_counts: Dict[str, int] = {}
    _seen_sample_names: set = set()
    samples: List[SampleInfo] = []
    idx = 0
    id_map_lines = [
        "internal_id\tsample_name\torig_token\ttopology\tsource_fasta\tinternal_fasta\torig_gb_path\n"
    ]
    for fa in in_files:
        try:
            records = list(SeqIO.parse(str(fa), "fasta"))
        except Exception as exc:
            logging.warning(f"[SPLIT] Cannot parse {fa}: {exc}")
            continue
        if not records:
            logging.warning(f"[SPLIT] No records in {fa}, skipping.")
            continue
        if (not if_split) and len(records) > 1:
            raise ValueError(f"FASTA has >1 record: {fa}. Use --if_split_fasta yes.")
        for rec in records:
            idx += 1
            desc = (rec.description or "").strip()
            orig_token = header_first_token(rec.id) if _is_topology_tag(desc) \
                         else header_first_token(desc)
            base_name = safe_name(orig_token)
            sample_counts[base_name] = sample_counts.get(base_name, 0) + 1
            cnt = sample_counts[base_name]
            cand_name = base_name if cnt == 1 else f"{base_name}_{cnt}"
            while cand_name in _seen_sample_names:
                sample_counts[base_name] += 1
                cand_name = f"{base_name}_{sample_counts[base_name]}"
            _seen_sample_names.add(cand_name)
            sample_name = cand_name
            topology    = infer_topology_from_description(desc)
            internal_id = make_internal_id(idx, id_maxlen)
            out_fa      = out_dir / f"{internal_id}.fasta"
            write_one_record_fasta(rec, out_fa, seq_id=internal_id, topology=topology)
            orig_gb_path = acc_to_orig_gb.get(safe_name(rec.id)) \
                        or acc_to_orig_gb.get(orig_token)
            si = SampleInfo(
                sample_name=sample_name, orig_token=orig_token,
                internal_id=internal_id, topology=topology,
                source_fasta=fa.resolve(), internal_fasta=out_fa.resolve(),
                orig_gb_path=orig_gb_path,
            )
            samples.append(si)
            id_map_lines.append(
                f"{internal_id}\t{sample_name}\t{orig_token}\t{topology}\t"
                f"{si.source_fasta}\t{si.internal_fasta}\t{orig_gb_path or 'NA'}\n"
            )
    id_map_path = out_dir / "id_map.tsv"
    id_map_path.write_text("".join(id_map_lines), encoding="utf-8")
    return samples, id_map_path


def parse_suffix_fq(s: str) -> Tuple[str, str]:
    parts = [p.strip() for p in (s or "").split(",") if p.strip()]
    if len(parts) != 2:
        raise ValueError("--suffix_fq needs two comma-separated suffixes")
    return parts[0], parts[1]


def find_fq_pair_for_sample(
    sample_name: str, infq: Path, fq_position: str, suffix1: str, suffix2: str,
) -> Tuple[Optional[Path], Optional[Path]]:
    sn = str(sample_name)
    if fq_position == "flat":
        f1 = infq / f"{sn}{suffix1}"
        f2 = infq / f"{sn}{suffix2}"
        return (f1.resolve() if f1.exists() else None), (f2.resolve() if f2.exists() else None)
    if fq_position == "sample_dir":
        d = infq / sn
        if not d.exists():
            return None, None
        f1 = d / f"{sn}{suffix1}"
        f2 = d / f"{sn}{suffix2}"
        if f1.exists() and f2.exists():
            return f1.resolve(), f2.resolve()
        c1 = sorted(d.glob(f"*{suffix1}"))
        c2 = sorted(d.glob(f"*{suffix2}"))
        return (c1[0].resolve() if c1 else None), (c2[0].resolve() if c2 else None)
    raise ValueError(f"Unknown --fq_position: {fq_position}")


def run_cmd_capture(
    cmd: List[str], cwd: Optional[Path], stdout_path: Path, stderr_path: Path,
) -> Tuple[int, float]:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    p  = subprocess.run(cmd, cwd=str(cwd) if cwd else None,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    dt = time.time() - t0
    if p.stdout:
        stdout_path.write_text(p.stdout, encoding="utf-8", errors="ignore")
    if p.stderr:
        stderr_path.write_text(p.stderr, encoding="utf-8", errors="ignore")
    return int(p.returncode), float(dt)


# ─────────────────────────────────────────────────────────────────────────────
# GenBank helpers
# ─────────────────────────────────────────────────────────────────────────────
def collect_genbank_files(src: Path) -> List[Path]:
    if src.is_file():
        return [src.resolve()]
    if src.is_dir():
        return sorted(
            p.resolve() for p in src.rglob("*")
            if p.is_file() and _is_gb_file(p)
        )
    return []


def merge_gb_files(gb_files: List[Path], merged_out: Path) -> int:
    merged_out.parent.mkdir(parents=True, exist_ok=True)
    all_recs: List["SeqRecord"] = []
    for gb in gb_files:
        try:
            all_recs.extend(list(SeqIO.parse(str(gb), "genbank")))
        except Exception as exc:
            logging.warning(f"[MERGE_SKIP] {gb}: {exc}")
    if all_recs:
        SeqIO.write(all_recs, str(merged_out), "genbank")
    return len(all_recs)


def find_genbank_like_files(root: Path, preferred_prefix: Optional[str] = None) -> List[Path]:
    out = [p for p in root.rglob("*")
           if p.is_file() and _is_gb_file(p)]
    if preferred_prefix:
        preferred = sorted([p for p in out if p.stem.startswith(preferred_prefix)],
                           key=lambda x: x.stat().st_size if x.exists() else 0, reverse=True)
        others    = sorted([p for p in out if not p.stem.startswith(preferred_prefix)],
                           key=lambda x: x.stat().st_size if x.exists() else 0, reverse=True)
        return preferred + others
    return sorted(out, key=lambda x: x.stat().st_size if x.exists() else 0, reverse=True)


def find_any_summary_txt(sample_dir: Path) -> Optional[Path]:
    cands = list(sample_dir.rglob("summary.txt"))
    if not cands:
        return None
    # BUG-V0.28-D: sort by actual directory depth (len(p.parts)) rather than
    # raw string length.  String length can be misleading when parent directory
    # names differ in length (e.g. a long flat name beating a deeply nested
    # short name), causing the wrong summary.txt to be selected.  MitoZ
    # always writes its summary.txt into a sub-subdirectory, so the deepest
    # file by parts count is the correct one; str(p) is used as a stable
    # tiebreak when two candidates share the same depth.
    cands.sort(key=lambda p: (len(p.parts), str(p)), reverse=True)
    return cands[0]


def parse_counts_from_summary_txt(summary_txt: Path) -> Optional[Tuple[int, int, int]]:
    txt = summary_txt.read_text(encoding="utf-8", errors="ignore")
    m1  = re.search(r"Protein\s+coding\s+genes\s+totally\s+found:\s*(\d+)", txt, re.I)
    m2  = re.search(r"tRNA\s+genes\s+totally\s+found:\s*(\d+)",            txt, re.I)
    m3  = re.search(r"rRNA\s+genes\s+totally\s+found:\s*(\d+)",            txt, re.I)
    if m1 and m2 and m3:
        return int(m1.group(1)), int(m2.group(1)), int(m3.group(1))
    return None


def counts_from_genbank(gb_path: Path) -> Tuple[int, int, int]:
    pcg = trna = rrna = 0
    for rec in SeqIO.parse(str(gb_path), "genbank"):
        for feat in rec.features:
            if feat.type == "CDS":    pcg  += 1
            elif feat.type == "tRNA": trna += 1
            elif feat.type == "rRNA": rrna += 1
    return pcg, trna, rrna


_LOCUS_NAME_MAXLEN = 16


def rewrite_genbank_record_ids(gb_in: Path, gb_out: Path, new_base_id: str) -> Tuple[int, int, int]:
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    if not recs:
        raise ValueError(f"No records in {gb_in}")
    pcg = trna = rrna = 0
    out_recs: List["SeqRecord"] = []
    for i, r in enumerate(recs, start=1):
        rid = new_base_id if i == 1 else f"{new_base_id}_{i}"
        # BUG-V0.28-C: r[:] triggers Biopython's __getitem__ which clones all
        # feature locations — wasteful because we immediately overwrite .id /
        # .name / .annotations anyway.  Build the new SeqRecord directly from
        # its constituent parts, mirroring the pattern in replace_organism_in_record
        # (BUG-NEW-E, v0.27).
        locus = rid[:_LOCUS_NAME_MAXLEN]
        if len(rid) > _LOCUS_NAME_MAXLEN:
            logging.info(
                f"[LOCUS_TRUNCATE] '{rid}' ({len(rid)} chars) truncated to '{locus}' "
                f"in LOCUS line (GenBank limit={_LOCUS_NAME_MAXLEN}). "
                f"ACCESSION/VERSION unaffected.")
        new_annotations = dict(r.annotations)
        new_annotations["accessions"] = [rid]
        new_annotations.pop("sequence_version", None)
        rr = SeqRecord(
            r.seq,
            id=rid,
            name=locus,
            description=r.description or "",
            dbxrefs=list(getattr(r, "dbxrefs", [])),
            annotations=new_annotations,
            features=list(r.features),
        )
        out_recs.append(rr)
        for feat in r.features:
            if feat.type == "CDS":    pcg  += 1
            elif feat.type == "tRNA": trna += 1
            elif feat.type == "rRNA": rrna += 1
    gb_out.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_recs, str(gb_out), "genbank")
    return pcg, trna, rrna


def gb_to_internal_fasta(gb_in: Path, fasta_out: Path, internal_id: str) -> None:
    """Convert a GenBank file to FASTA for Phase-2 re-annotation input.

    All records are written (not just recs[0]); secondary contigs get
    suffixes _2, _3, … appended to internal_id (v0.24 / Hole-1).
    """
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    if not recs:
        raise ValueError(f"No records in {gb_in}")
    out_recs = []
    for i, rec in enumerate(recs):
        topo   = rec.annotations.get("topology", "circular")
        seq_id = internal_id if i == 0 else f"{internal_id}_{i + 1}"
        out_recs.append(SeqRecord(rec.seq, id=seq_id, name=seq_id,
                                  description=f"topology={topo}"))
    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    with fasta_out.open("w", encoding="utf-8") as fh:
        SeqIO.write(out_recs, fh, "fasta")


# ─────────────────────────────────────────────────────────────────────────────
# GenBank reorientation
# ─────────────────────────────────────────────────────────────────────────────
def normalize_gene_key(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (name or "").strip().lower())


def gene_aliases(gene: str) -> List[str]:
    g = normalize_gene_key(gene)
    if g in ("cox1", "coi", "coxi"):
        return ["cox1", "coi", "coxi",
                "cytochromecoxidasesubuniti", "cytochromecoxidasesubunit1"]
    return [g]


def feature_gene_names(feature: "SeqFeature") -> List[str]:
    names: List[str] = []
    for k in ("gene", "product", "label"):
        for v in feature.qualifiers.get(k, []):
            names.append(str(v))
    return names


def find_gene_feature(record: "SeqRecord", target_gene: str) -> Optional["SeqFeature"]:
    aliases = set(gene_aliases(target_gene))
    for feat in record.features:
        if feat.type not in {"gene", "CDS", "tRNA", "rRNA"}:
            continue
        for nm in feature_gene_names(feat):
            if normalize_gene_key(nm) in aliases:
                return feat
    return None


def _rc_location(loc, seq_len: int):
    if isinstance(loc, CompoundLocation):
        return CompoundLocation(
            [_rc_location(p, seq_len) for p in loc.parts][::-1], operator=loc.operator)
    start, end = int(loc.start), int(loc.end)
    new_strand  = None if loc.strand is None else (-loc.strand)
    return FeatureLocation(seq_len - end, seq_len - start, strand=new_strand)


def reverse_complement_record(record: "SeqRecord") -> "SeqRecord":
    seq_len = len(record.seq)
    new_rec = SeqRecord(record.seq.reverse_complement(), id=record.id,
                        name=record.name, description=record.description)
    new_rec.annotations = dict(record.annotations)
    new_rec.dbxrefs     = list(getattr(record, "dbxrefs", []))
    new_features: List["SeqFeature"] = []
    for feat in record.features:
        try:
            new_loc = _rc_location(feat.location, seq_len)
            new_features.append(SeqFeature(location=new_loc, type=feat.type,
                                            id=getattr(feat, "id", "<unknown>"),
                                            qualifiers=dict(feat.qualifiers)))
        except Exception:
            new_features.append(feat)
    new_rec.features = new_features
    return new_rec


def shift_location(location, shift: int, seq_len: int):
    """Shift all feature coordinates by -shift (mod seq_len).

    BUG-FIX (v0.20): Flatten nested CompoundLocations.
    BUG-FIX (v0.21a): Ghost zero-length feature at cut-point.
    BUG-FIX (v0.21b): Genome-spanning source feature shattering.
    BUG-FIX (v0.22): Use <= for zero-length point features.
    """
    if isinstance(location, CompoundLocation):
        flat_parts: List = []
        for p in location.parts:
            shifted = shift_location(p, shift, seq_len)
            if isinstance(shifted, CompoundLocation):
                flat_parts.extend(shifted.parts)
            else:
                flat_parts.append(shifted)
        return CompoundLocation(flat_parts, operator=location.operator)

    orig_start = int(location.start)
    orig_end   = int(location.end)
    strand     = location.strand
    length     = orig_end - orig_start

    if length == seq_len:
        return FeatureLocation(0, seq_len, strand=strand)

    start = (orig_start - shift) % seq_len
    end   = (orig_end   - shift) % seq_len

    if end == 0 and length > 0:
        end = seq_len

    if start <= end:
        return FeatureLocation(start, end, strand=strand)

    return CompoundLocation([FeatureLocation(start, seq_len, strand=strand),
                              FeatureLocation(0, end, strand=strand)])


def rotate_record(record: "SeqRecord", start_pos: int) -> "SeqRecord":
    seq_len = len(record.seq)
    topology = record.annotations.get("topology", "circular")
    if str(topology).lower() == "linear" and start_pos != 0:
        logging.warning(
            f"[LINEAR_REORIENT_SKIP] record={record.id}: target start={start_pos}; "
            "linear topology cannot be circularly rotated, sequence left unchanged.")
        # BUG-V0.29-M: never circularly rotate a linear contig.  The previous
        # implementation preserved topology but still moved the sequence tail
        # before the head, which corrupts incomplete linear contigs.
        # new_seq = record.seq[start_pos:] + record.seq[:start_pos]
        return SeqRecord(
            record.seq,
            id=record.id, name=record.name, description=record.description,
            dbxrefs=list(getattr(record, "dbxrefs", [])),
            annotations=dict(record.annotations),
            features=list(record.features),
        )
    new_seq = record.seq[start_pos:] + record.seq[:start_pos]
    new_rec = SeqRecord(
        Seq(str(new_seq)),
        id=record.id, name=record.name, description=record.description,
        dbxrefs=list(getattr(record, "dbxrefs", [])),
        # BUG-V0.29-B: the old line forced every reoriented record to
        # circular, which changed true linear contigs into circular records.
        # annotations={**record.annotations, "topology": "circular"},
        annotations={**record.annotations, "topology": topology},
    )
    new_features: List["SeqFeature"] = []
    for feat in record.features:
        try:
            new_loc = shift_location(feat.location, start_pos, seq_len)
            new_features.append(SeqFeature(location=new_loc, type=feat.type,
                                            id=getattr(feat, "id", "<unknown>"),
                                            qualifiers=dict(feat.qualifiers)))
        except Exception:
            new_features.append(feat)
    new_features.sort(key=lambda f: feature_start_pos(f))
    new_rec.features = new_features
    return new_rec


def feature_start_pos(feat: "SeqFeature") -> int:
    loc = feat.location
    if isinstance(loc, CompoundLocation) and loc.parts:
        return int(loc.parts[0].start)
    return int(loc.start)


def reorient_one_gb(
    gb_in: Path, gb_out: Path, gb_no_gene: Path,
    which_gene_first: str, force_forward: bool,
) -> Tuple[int, int, List[str]]:
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    if not recs:
        raise ValueError(f"No GenBank records: {gb_in}")
    ok: List["SeqRecord"] = []
    no_gene: List["SeqRecord"] = []
    no_gene_ids: List[str] = []
    for rec in recs:
        try:
            feat = find_gene_feature(rec, which_gene_first)
            if feat is None:
                no_gene.append(rec); no_gene_ids.append(rec.id); continue
            work_rec = rec
            if force_forward and feat.location.strand == -1:
                work_rec = reverse_complement_record(rec)
                feat = find_gene_feature(work_rec, which_gene_first)
                if feat is None:
                    no_gene.append(rec); no_gene_ids.append(rec.id); continue
            ok.append(rotate_record(work_rec, feature_start_pos(feat)))
        except Exception as exc:
            logging.warning(f"[REORIENT_REC_FAIL] record={rec.id}: {exc}")
            no_gene.append(rec); no_gene_ids.append(rec.id)
    gb_out.parent.mkdir(parents=True, exist_ok=True)
    gb_no_gene.parent.mkdir(parents=True, exist_ok=True)
    if ok:
        SeqIO.write(ok, str(gb_out), "genbank")
    if no_gene:
        SeqIO.write(no_gene, str(gb_no_gene), "genbank")
    return len(ok), len(no_gene), no_gene_ids


# ─────────────────────────────────────────────────────────────────────────────
# Metadata transfer engine
# ─────────────────────────────────────────────────────────────────────────────
_REFSEQ_PREFIX_PAT = re.compile(r"^([A-Z]{2}_)")
# BUG-V0.29-A: Only accession-shaped IDs should have a trailing version
# stripped.  Local sample names such as B64_2 are distinct samples, not
# accession versions.
_NON_REFSEQ_VERSION_PAT = re.compile(r"^[A-Z]{1,4}\d{5,}[._]\d+$", re.I)
# Pre-compiled regex for LOCUS bp-count substitution (OPT-1 v0.27)
_LOCUS_BPCOUNT_RE  = re.compile(r"(LOCUS\s+\S+\s+)\d+(\s+[a-zA-Z]+)")


def _clean_acc(acc: Optional[str]) -> Optional[str]:
    if not acc:
        return None
    return acc.strip().rstrip(";").strip().replace(" ", "") or None


def _normalize_acc(acc: Optional[str]) -> Optional[str]:
    s = _clean_acc(acc)
    if not s:
        return None
    m = _REFSEQ_PREFIX_PAT.match(s)
    if m:
        prefix = m.group(1)
        body   = re.sub(r"[._]\d+$", "", s[len(prefix):])
        s      = prefix + body
    else:
        # BUG-V0.29-A: the old broad substitution below collapsed local IDs
        # such as B64_2 -> B64.  Keep it as a comment for auditability.
        # s = re.sub(r"[._]\d+$", "", s)
        if _NON_REFSEQ_VERSION_PAT.match(s):
            s = re.sub(r"[._]\d+$", "", s)
    return s.rstrip(";").strip() or None


@functools.lru_cache(maxsize=4096)
def _acc_variants(acc: Optional[str]) -> Tuple[str, ...]:
    """Return all meaningful accession variants for lookup matching.

    OPT-2 (v0.27): cached with lru_cache — repeated calls with the same
    accession string (common during index building) are free after the
    first computation.

    BUG-V0.28-A: return type changed from List[str] to Tuple[str, ...] so
    the cached object is immutable.  The previous List was safe only by
    convention ("callers must not mutate the result"), which is fragile;
    any accidental append/clear/sort on the returned list would silently
    corrupt the cache.  All call-sites iterate over the result, so the
    change is fully backward-compatible.
    """
    clean = _clean_acc(acc)
    if not clean:
        return ()
    norm = _normalize_acc(clean)
    seen: List[str] = []

    def _add(x: Optional[str]) -> None:
        if x and x not in seen:
            seen.append(x)

    _add(clean)
    _add(norm)
    m = re.search(r"[._](\d+)$", clean)
    # BUG-V0.29-A: only add version-style alternatives when normalisation
    # actually removed a version suffix.  Otherwise B64_2 would gain B64 as
    # an alias and collide with the real B64 sample.
    # if m and norm:
    if m and norm and norm != clean:
        ver = m.group(1)
        _add(f"{norm}.{ver}"); _add(f"{norm}_{ver}")
    return tuple(seen)


def _split_gb_records(text: str) -> List[str]:
    """Split a multi-record GenBank text into individual record strings.

    BUG-NEW-C (v0.27): empty / whitespace-only pseudo-records that result
    from trailing newlines after the final ``//`` terminator are now
    filtered out, preventing phantom entries in replace_metadata_in_text
    stats counters.
    """
    records: List[str] = []
    buf: List[str] = []
    for line in text.splitlines(keepends=True):
        buf.append(line)
        if line.strip() == "//":
            chunk = "".join(buf)
            if chunk.strip():          # skip empty / whitespace-only chunks
                records.append(chunk)
            buf = []
    if buf:
        chunk = "".join(buf)
        if chunk.strip():
            records.append(chunk)
    return records


def _split_meta_body(record: str) -> Tuple[str, str]:
    lines = record.splitlines(keepends=True)
    meta: List[str] = []
    body: List[str] = []
    in_body = False
    for line in lines:
        if not in_body and (line.startswith("FEATURES") or line.startswith("ORIGIN")):
            in_body = True
        (body if in_body else meta).append(line)
    return "".join(meta), "".join(body)


def _get_record_accession(record: str) -> Optional[str]:
    """Extract accession from GB text — VERSION first, then ACCESSION, then LOCUS.

    BUG-FIX (v0.22): Biopython truncates LOCUS name to 16 chars; scanning
    VERSION / ACCESSION first avoids silently returning a truncated key.
    """
    for prefix in ("VERSION", "ACCESSION", "LOCUS"):
        for line in record.splitlines():
            if line.startswith(prefix):
                parts = line.split()
                if len(parts) >= 2:
                    return _clean_acc(parts[1])
    return None


def _meta_is_meaningful(meta: str) -> bool:
    if not meta.strip():
        return False
    lines = [ln for ln in meta.splitlines() if ln.strip()]
    if not lines or not lines[0].startswith("LOCUS"):
        return False
    useful = ("DEFINITION", "ACCESSION", "VERSION", "KEYWORDS",
              "SOURCE", "REFERENCE", "COMMENT", "DBLINK", "ORGANISM")
    return any(any(ln.startswith(p) for p in useful) for ln in lines[1:])


def build_metadata_index(src_text: str) -> Tuple[Dict, Dict, Dict]:
    exact_dict: Dict[str, str] = {}
    clean_dict: Dict[str, str] = {}
    norm_dict:  Dict[str, Tuple[str, str]] = {}
    for rec in _split_gb_records(src_text):
        raw_acc = _get_record_accession(rec)
        if not raw_acc:
            continue
        meta, _ = _split_meta_body(rec)
        if not _meta_is_meaningful(meta):
            continue
        for v in _acc_variants(raw_acc):
            exact_dict[v] = meta
        cl = _clean_acc(raw_acc)
        if cl:
            clean_dict.setdefault(cl, meta)
        norm = _normalize_acc(raw_acc)
        if norm:
            norm_dict.setdefault(norm, (raw_acc, meta))
    return exact_dict, clean_dict, norm_dict


def find_metadata_from_index(acc: str, exact_dict, clean_dict, norm_dict,
                              ) -> Tuple[Optional[str], Optional[str]]:
    for v in _acc_variants(acc):
        if v in exact_dict:
            mtype = "exact" if v == acc else f"clean '{acc}'->'{v}'"
            return exact_dict[v], mtype
    cl = _clean_acc(acc)
    if cl and cl in clean_dict:
        return clean_dict[cl], f"clean '{acc}'->'{cl}'"
    norm = _normalize_acc(acc)
    if norm and norm in norm_dict:
        orig, meta = norm_dict[norm]
        return meta, f"compat '{orig}'<->'{acc}'"
    return None, None


def _similarity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    return difflib.SequenceMatcher(None, a.upper(), b.upper()).ratio()


def diagnose_unmatched(acc: str, exact_dict: Dict) -> Tuple[List[str], List[Tuple[int, str]]]:
    hints: List[str] = []
    candidates: List[Tuple[int, str]] = []
    acc_cl   = _clean_acc(acc)     or acc
    acc_norm = _normalize_acc(acc) or acc
    acc_up   = acc_cl.upper()
    seen: set = set()
    for k in exact_dict:
        ck = _clean_acc(k)
        if not ck or ck in seen:
            continue
        seen.add(ck)
        orig_norm = _normalize_acc(ck) or ck
        score = 0; reason: Optional[str] = None
        if ck.upper() == acc_up and ck != acc_cl:
            score = 90; reason = f"[case mismatch] orig='{k}' new='{acc}'"
        elif orig_norm == acc_norm and ck != acc_cl:
            score = 88; reason = f"[version/style mismatch] orig='{k}' new='{acc}'"
        else:
            sim = _similarity(acc_norm, orig_norm)
            if sim > 0.78:
                score = int(sim * 55)
        if score > 0:
            if reason:
                hints.append(reason)
            candidates.append((score, k))
    if not hints:
        hints.append(f"[no diagnosis] no similar accession found for '{acc}'")
    candidates.sort(key=lambda x: -x[0])
    return hints[:5], candidates[:5]


def interactive_resolve(acc: str, candidates: List[Tuple[int, str]],
                        exact_dict: Dict, auto_accept: bool = False,
                        ) -> Tuple[Optional[str], Optional[str]]:
    if not candidates:
        return None, None
    best_score, best_acc = candidates[0]
    if auto_accept:
        if best_score >= 60:
            logging.info(f"[AUTO-ACCEPT] '{acc}' -> '{best_acc}' (score={best_score})")
            return exact_dict.get(best_acc), best_acc
        logging.warning(f"[AUTO-SKIP] '{acc}' top score {best_score} < 60")
        return None, None
    if not sys.stdin.isatty():
        logging.warning(
            f"[NON-INTERACTIVE] stdin is not a TTY; skipping interactive resolution "
            f"for '{acc}'. Top candidate: '{best_acc}' (score={best_score}). "
            f"Use --auto_accept_metadata yes to resolve automatically.")
        return None, None
    print(f"\n⚠ '{acc}' not matched. Candidates:")
    for i, (score, orig_acc) in enumerate(candidates, start=1):
        print(f"  [{i}] {orig_acc:<30} score={score}")
    print("  [s] skip   [q] quit")
    while True:
        try:
            choice = input("choose: ").strip().lower()
        except (EOFError, KeyboardInterrupt):
            print(); return None, None
        if choice == "q":
            sys.exit(0)
        if choice == "s":
            return None, None
        if choice.isdigit():
            idx = int(choice) - 1
            if 0 <= idx < len(candidates):
                chosen = candidates[idx][1]
                return exact_dict.get(chosen), chosen
        print("  invalid, retry.")


def replace_metadata_in_text(
    new_gb_text: str, exact_dict, clean_dict, norm_dict,
    smart: bool = False, auto_accept: bool = False,
) -> Tuple[str, Dict]:
    """Replace GB header sections in *new_gb_text* using the metadata index.

    OPT-1 (v0.27): the LOCUS bp-count substitution regex is pre-compiled
    at module level (_LOCUS_BPCOUNT_RE) rather than recompiled per record.
    """
    stats: Dict = {
        "total_new": 0, "no_accession": 0,
        "exact": 0, "compat": 0, "smart_resolved": 0, "not_found": 0,
        "exact_list": [], "compat_list": [], "smart_list": [], "not_found_list": [],
    }
    out_records: List[str] = []
    for rec in _split_gb_records(new_gb_text):
        stats["total_new"] += 1
        acc = _get_record_accession(rec)
        if not acc:
            stats["no_accession"] += 1
            out_records.append(rec); continue
        meta, match_type = find_metadata_from_index(acc, exact_dict, clean_dict, norm_dict)
        if meta is not None:
            _, body = _split_meta_body(rec)
            # Fix the bp count in the source metadata's LOCUS line to match
            # the actual (re-annotated) sequence length (BUG-3 v0.20).
            new_locus_m = re.search(r"LOCUS\s+\S+\s+(\d+)\s+[a-zA-Z]+", rec)
            if new_locus_m:
                real_len = new_locus_m.group(1)
                meta = _LOCUS_BPCOUNT_RE.sub(rf"\g<1>{real_len}\g<2>", meta, count=1)
            out_records.append(meta + body if body else rec)
            if match_type == "exact":
                stats["exact"] += 1; stats["exact_list"].append(acc)
            else:
                stats["compat"] += 1; stats["compat_list"].append((acc, match_type))
            continue
        hints, candidates = diagnose_unmatched(acc, exact_dict)
        if smart and candidates:
            resolved_meta, chosen_acc = interactive_resolve(
                acc, candidates, exact_dict, auto_accept=auto_accept)
            if resolved_meta is not None:
                _, body = _split_meta_body(rec)
                out_records.append(resolved_meta + body if body else rec)
                stats["smart_resolved"] += 1
                stats["smart_list"].append((acc, chosen_acc))
                continue
        stats["not_found"] += 1
        stats["not_found_list"].append((acc, hints))
        out_records.append(rec)
    return "".join(out_records), stats


def transfer_gb_metadata_file(
    metadata_src_gb: Path, new_body_gb: Path, output_gb: Path,
    smart: bool = False, auto_accept: bool = False,
) -> Tuple[bool, int, int, List[str], str, List[str]]:
    if not metadata_src_gb.exists():
        return False, 0, 0, [], f"metadata source not found: {metadata_src_gb}", []
    if not new_body_gb.exists():
        return False, 0, 0, [], f"body source not found: {new_body_gb}", []
    src_text = metadata_src_gb.read_text(encoding="utf-8", errors="replace")
    new_text = new_body_gb.read_text(encoding="utf-8", errors="replace")
    exact_dict, clean_dict, norm_dict = build_metadata_index(src_text)
    if not exact_dict and not norm_dict:
        return False, 0, 0, [], "metadata source has no usable metadata", []
    replaced_text, stats = replace_metadata_in_text(
        new_text, exact_dict, clean_dict, norm_dict, smart=smart, auto_accept=auto_accept)
    output_gb.parent.mkdir(parents=True, exist_ok=True)
    output_gb.write_text(replaced_text, encoding="utf-8")
    n_replaced  = stats["exact"] + stats["compat"] + stats["smart_resolved"]
    unmatched   = [x[0] for x in stats["not_found_list"]]
    compat_keys = [x[0] for x in stats["compat_list"]]
    return True, n_replaced, stats["not_found"], unmatched, "ok", compat_keys


# ─────────────────────────────────────────────────────────────────────────────
# Organism / SOURCE replacement helpers  (v0.25 / v0.26 / v0.27)
# ─────────────────────────────────────────────────────────────────────────────
def build_organism_index(gb_files: List[Path]) -> Dict[str, Tuple[str, str]]:
    """Build accession → (organism_name, source_line_text) lookup.

    BUG-NEW-B (v0.27): removed the redundant ``index.setdefault(norm,
    entry)`` executed after the _acc_variants loop.  _acc_variants already
    registers the normalised form as one of its variants; the extra call
    was dead code and could inadvertently mask a later (higher-priority)
    entry under the same norm key.
    """
    index: Dict[str, Tuple[str, str]] = {}
    for gb_path in gb_files:
        try:
            records = list(SeqIO.parse(str(gb_path), "genbank"))
        except Exception as exc:
            logging.warning(f"[ORG_INDEX] Cannot parse {gb_path}: {exc}")
            continue
        for rec in records:
            acc = gb_record_accession(rec)

            organism: str = rec.annotations.get("organism", "").strip()
            if not organism:
                for feat in rec.features:
                    if feat.type == "source":
                        vals = feat.qualifiers.get("organism", [])
                        if vals:
                            organism = str(vals[0]).strip()
                        break
            if not organism:
                organism = rec.annotations.get("source", "").strip()

            source_text: str = rec.annotations.get("source", "").strip()

            if not organism and not source_text:
                continue

            entry: Tuple[str, str] = (organism, source_text)

            # Register under all accession variants for robust matching.
            # _acc_variants already includes the normalised form, so no
            # additional setdefault(norm, entry) call is needed here.
            for variant in _acc_variants(acc):
                # BUG-V0.29-K: old setdefault silently kept the first
                # organism/source when two metadata records normalised to the
                # same accession variant.  Preserve first-wins behaviour but
                # make the conflict visible.
                # index.setdefault(variant, entry)
                if variant in index and index[variant] != entry:
                    logging.warning(
                        f"[ORG_INDEX_CONFLICT] {variant}: keeping first "
                        f"{index[variant]!r}, ignoring {entry!r} from {gb_path.name}")
                else:
                    index.setdefault(variant, entry)

    return index


def lookup_organism(
    acc: str, org_index: Dict[str, Tuple[str, str]]
) -> Optional[Tuple[str, str]]:
    """Return (organism, source_text) for *acc* using exact → normalised fallback."""
    for v in _acc_variants(acc):
        if v in org_index:
            return org_index[v]
    norm = _normalize_acc(acc)
    if norm and norm in org_index:
        return org_index[norm]
    return None


def replace_organism_in_record(
    rec: "SeqRecord",
    new_organism: str,
    new_source: Optional[str] = None,
) -> "SeqRecord":
    """Return a copy of *rec* with the organism name replaced throughout.

    BUG-NEW-E (v0.27): eliminated the wasteful rec[:] slice-copy.  The
    old code created a full Biopython slice (copying all features with
    original locations) only to immediately overwrite new_rec.features
    with the reprocessed list.  The SeqRecord is now built directly from
    its constituent parts.

    Updates:
      • rec.annotations["organism"]  — always set when new_organism is non-empty.
      • rec.annotations["source"]    — set only when new_source is explicitly
                                       provided and non-empty (BUG-11 v0.26).
      • /organism= qualifier on every ``source`` feature — unconditionally
                                       set (BUG-A v0.26).
    """
    new_annotations = dict(rec.annotations)
    if new_organism:
        new_annotations["organism"] = new_organism
    if new_source:
        new_annotations["source"] = new_source

    new_features: List["SeqFeature"] = []
    for feat in rec.features:
        new_quals = dict(feat.qualifiers)
        if feat.type == "source" and new_organism:
            new_quals["organism"] = [new_organism]
        new_features.append(
            SeqFeature(
                location=feat.location,
                type=feat.type,
                id=getattr(feat, "id", "<unknown>"),
                qualifiers=new_quals,
            )
        )

    new_rec = SeqRecord(
        rec.seq,
        id=rec.id, name=rec.name, description=rec.description,
        dbxrefs=list(getattr(rec, "dbxrefs", [])),
        annotations=new_annotations,
        features=new_features,
    )
    return new_rec


def replace_organism_in_gb_file(
    gb_in:     Path,
    gb_out:    Path,
    org_index: Dict[str, Tuple[str, str]],
) -> Tuple[int, int, int]:
    """Replace /organism= (and optionally SOURCE) in every record of *gb_in*."""
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    if not recs:
        raise ValueError(f"No records in {gb_in}")
    out_recs: List["SeqRecord"] = []
    n_replaced = n_not_found = 0
    for rec in recs:
        acc    = gb_record_accession(rec)
        result = lookup_organism(acc, org_index)
        if result is not None:
            organism, source_text = result
            out_recs.append(
                replace_organism_in_record(rec, organism,
                                           new_source=source_text if source_text else None)
            )
            n_replaced += 1
            logging.info(
                f"[REPLACE_ORG] {acc}: /organism= → '{organism}'"
                + (f"  SOURCE → '{source_text}'" if source_text else "")
            )
        else:
            out_recs.append(rec)
            n_not_found += 1
            logging.warning(
                f"[REPLACE_ORG_MISS] {acc}: no entry in organism index — unchanged.")
    gb_out.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_recs, str(gb_out), "genbank")
    return n_replaced, n_not_found, len(recs)


def _sync_organism_qualifier(gb_in: Path, gb_out: Path) -> int:
    """Sync /organism= in FEATURES source to match rec.annotations["organism"].

    BUG-NEW-E (v0.27): same rec[:] waste removed as in
    replace_organism_in_record — SeqRecord built directly.
    BUG-B (v0.26): qualifier is added even when previously absent.
    """
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    out_recs: List["SeqRecord"] = []
    n_fixed = 0
    for rec in recs:
        organism = rec.annotations.get("organism", "").strip()
        if not organism:
            out_recs.append(rec)
            continue
        new_features: List["SeqFeature"] = []
        changed = False
        for feat in rec.features:
            new_quals = dict(feat.qualifiers)
            if feat.type == "source":
                current = str(new_quals["organism"][0]).strip() \
                          if new_quals.get("organism") else ""
                if current != organism:
                    new_quals["organism"] = [organism]
                    changed = True
            new_features.append(
                SeqFeature(
                    location=feat.location,
                    type=feat.type,
                    id=getattr(feat, "id", "<unknown>"),
                    qualifiers=new_quals,
                )
            )
        new_rec = SeqRecord(
            rec.seq,
            id=rec.id, name=rec.name, description=rec.description,
            dbxrefs=list(getattr(rec, "dbxrefs", [])),
            annotations=dict(rec.annotations),
            features=new_features,
        )
        out_recs.append(new_rec)
        if changed:
            n_fixed += 1
    gb_out.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_recs, str(gb_out), "genbank")
    return n_fixed


def gb_file_has_meaningful_metadata(gb_path: Path) -> bool:
    if not gb_path.is_file():
        return False
    try:
        txt = gb_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return False
    return any(_meta_is_meaningful(_split_meta_body(r)[0]) for r in _split_gb_records(txt))


# ─────────────────────────────────────────────────────────────────────────────
# Environment preflight
# ─────────────────────────────────────────────────────────────────────────────
def check_environment(args: argparse.Namespace) -> List[str]:
    problems: List[str] = []
    if SeqIO is None:
        raise RuntimeError("Biopython not installed. Run: pip install biopython")
    try:
        import Bio
        logging.info(f"[ENV] Biopython {getattr(Bio, '__version__', '?')} ✔")
    except Exception:
        pass
    major, minor = sys.version_info.major, sys.version_info.minor
    if (major, minor) < (3, 8):
        raise RuntimeError(f"Python 3.8+ required, found {major}.{minor}")
    logging.info(f"[ENV] Python {major}.{minor} ✔")

    if getattr(args, "cmd", None) in ("annotate", "run"):
        conda_env = getattr(args, "conda_env_mitoz", None)
        tool_path = getattr(args, "mitoz_path", None)
        conda_exe = getattr(args, "conda_exe", "conda")
        try:
            test_cmd = build_tool_cmd(
                tool_label="MitoZ", exe_name=MITOZ_EXE_NAME, base_args=["--version"],
                conda_env=conda_env, tool_path=tool_path, conda_exe=conda_exe)
        except FileNotFoundError as exc:
            raise RuntimeError(f"[ENV] MitoZ not reachable: {exc}") from exc
        logging.info(f"[ENV] Testing MitoZ: {' '.join(test_cmd)}")
        try:
            result = subprocess.run(test_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    text=True, timeout=30)
            combined = (result.stdout or "") + (result.stderr or "")
            if result.returncode == 0 or "mitoz" in combined.lower():
                ver_line = next((ln.strip() for ln in combined.splitlines() if ln.strip()), "?")
                logging.info(f"[ENV] MitoZ reachable ✔  ({ver_line[:80]})")
            else:
                problems.append(f"[ENV] MitoZ test rc={result.returncode}. Output: {combined[:200]}")
        except subprocess.TimeoutExpired:
            problems.append("[ENV] MitoZ test timed out after 30s")
        except Exception as exc:
            raise RuntimeError(f"[ENV] MitoZ execution error: {exc}") from exc

        if conda_env:
            conda_bin = which_or_none(conda_exe)
            if conda_bin:
                try:
                    r = subprocess.run([conda_bin, "env", "list"],
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       text=True, timeout=15)
                    env_names = [ln.split()[0] for ln in r.stdout.splitlines()
                                 if ln.strip() and not ln.startswith("#")]
                    if conda_env not in env_names:
                        problems.append(
                            f"[ENV] Conda env '{conda_env}' not found. "
                            f"Available: {', '.join(env_names[:10])}")
                    else:
                        logging.info(f"[ENV] Conda env '{conda_env}' found ✔")
                except Exception as exc:
                    problems.append(f"[ENV] Could not list conda envs: {exc}")
            else:
                problems.append(
                    f"[ENV] '{conda_exe}' not in PATH but --conda_env_mitoz='{conda_env}' set.")

    smk_bin = getattr(args, "snakemake_bin", None)
    backend = getattr(args, "scheduler_backend", "threads")
    if backend in ("snakemake", "auto"):
        if snakemake_available(smk_bin):
            smk_exe = smk_bin or which_or_none("snakemake") or "snakemake"
            try:
                sv = subprocess.run([smk_exe, "--version"], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, text=True, timeout=10)
                smk_ver = (sv.stdout or sv.stderr or "").strip().splitlines()[0] \
                          if sv.returncode == 0 else "?"
                logging.info(f"[ENV] Snakemake {smk_ver} ✔")
            except Exception:
                logging.info("[ENV] Snakemake available (version check failed)")
        else:
            if backend == "snakemake":
                raise RuntimeError("[ENV] --scheduler_backend snakemake requested but snakemake not found.")
            logging.info("[ENV] Snakemake not found, will use threads backend.")

    if hasattr(args, "input") and args.input:
        in_p = Path(args.input)
        if not in_p.exists():
            raise RuntimeError(f"[ENV] Input path does not exist: {in_p}")
        logging.info(f"[ENV] Input path exists ✔  ({in_p})")

    if hasattr(args, "infq") and args.infq:
        fq_p = Path(args.infq)
        if not fq_p.exists():
            problems.append(f"[ENV] --infq path does not exist: {fq_p}")
        else:
            logging.info(f"[ENV] FASTQ dir exists ✔  ({fq_p})")

    return problems


def preflight(args: argparse.Namespace) -> None:
    problems = check_environment(args)
    if problems:
        for p in problems:
            logging.warning(p)
        logging.warning(f"[ENV] {len(problems)} warning(s) above — pipeline continues but may fail.")


# ─────────────────────────────────────────────────────────────────────────────
# MitoZ annotation worker
# ─────────────────────────────────────────────────────────────────────────────
def mitoz_annotate_one(
    si: SampleInfo, fasta_file: Path, out_root: Path, phase_dirname: str,
    args: argparse.Namespace, fq1: Optional[Path], fq2: Optional[Path],
    total: int = 0, done_counter: Optional[List] = None,
) -> Tuple[str, bool, Path, Optional[Path], Optional[Path], Tuple[int, int, int]]:
    internal_id = si.internal_id
    workdir     = out_root / phase_dirname / internal_id
    workdir.mkdir(parents=True, exist_ok=True)
    log_name = f"sample.{phase_dirname}.{internal_id}"
    lg = get_sample_logger(log_name, workdir / "sample.log")
    lg.info(f"[SAMPLE] internal_id={internal_id} sample_name={si.sample_name}")
    lg.info(f"[INPUT_FASTA] {fasta_file}")

    base_args: List[str] = [
        "annotate",
        "--workdir",       str(workdir),
        "--outprefix",     internal_id,
        "--thread_number", str(int(args.threads)),
        "--clade",         str(args.clade),
        "--genetic_code",  str(int(args.genetic_code)),
        "--fastafiles",    str(fasta_file),
    ]
    if args.species_name:
        base_args += ["--species_name", str(args.species_name)]
    if fq1 and fq2:
        base_args += ["--fq1", str(fq1), "--fq2", str(fq2)]
        lg.info(f"[FQ] fq1={fq1}  fq2={fq2}")
    for x in args.mitoz_extra or []:
        base_args += shlex.split(str(x))

    cmd = build_tool_cmd(
        tool_label="MitoZ", exe_name=MITOZ_EXE_NAME, base_args=base_args,
        conda_env=args.conda_env_mitoz, tool_path=args.mitoz_path,
        conda_exe=args.conda_exe)
    lg.info("[CMD] " + " ".join(cmd))

    logs_dir    = workdir / "logs"
    stdout_path = logs_dir / "mitoz_annotate.stdout.log"
    stderr_path = logs_dir / "mitoz_annotate.stderr.log"
    rc, dt      = run_cmd_capture(cmd, cwd=workdir, stdout_path=stdout_path, stderr_path=stderr_path)
    if rc != 0:
        (logs_dir / "mitoz_annotate.returncode.txt").write_text(str(rc) + "\n", encoding="utf-8")
    (logs_dir / "mitoz_annotate.runtime_sec.txt").write_text(f"{dt:.2f}\n", encoding="utf-8")

    if done_counter is not None:
        lock: threading.Lock = done_counter[1]
        with lock:
            done_counter[0] += 1
            done_n = done_counter[0]
        logging.info(f"[PROGRESS] {done_n}/{total} '{si.sample_name}' done (rc={rc}, {dt:.0f}s)")

    if rc != 0:
        lg.error(f"[FAIL] MitoZ annotate returned {rc}")
        lg.error("---- stderr tail ----\n" +
                 tail_text(stderr_path.read_text(encoding="utf-8", errors="ignore")
                           if stderr_path.exists() else ""))
        close_sample_logger(log_name)
        return internal_id, False, workdir, None, None, (0, 0, 0)

    gb_files    = find_genbank_like_files(workdir, preferred_prefix=internal_id)
    best_gb     = gb_files[0] if gb_files else None
    summary_txt = find_any_summary_txt(workdir)

    counts: Optional[Tuple[int, int, int]] = None
    if summary_txt:
        counts = parse_counts_from_summary_txt(summary_txt)
    if counts is None and best_gb:
        try:
            counts = counts_from_genbank(best_gb)
        except Exception:
            counts = (0, 0, 0)
    counts = counts or (0, 0, 0)

    close_sample_logger(log_name)
    return internal_id, True, workdir, best_gb, summary_txt, counts


# ─────────────────────────────────────────────────────────────────────────────
# Snakemake scheduler
# ─────────────────────────────────────────────────────────────────────────────
def snakemake_available(snakemake_bin: Optional[str]) -> bool:
    if snakemake_bin:
        return Path(snakemake_bin).exists() or which_or_none(snakemake_bin) is not None
    return which_or_none("snakemake") is not None


@dataclass
class AnnotTask:
    internal_id:   str
    sample_name:   str
    fasta_file:    str
    phase_dirname: str
    out_root:      str
    fq1:           Optional[str]
    fq2:           Optional[str]


def render_snakemake_worker_script(out_dir: Path) -> Path:
    """Generate the self-contained MitoZ annotation worker script.

    NOTE (v0.27): This worker handles MitoZ annotation only.
    finalize_genbank_output (ID-rewrite, reorientation, metadata overlay,
    /organism= sync) is always executed in the main process after all
    Snakemake jobs complete, so all pipeline bug-fixes are applied
    regardless of the scheduler backend.
    """
    fp   = out_dir / "_mitoz_worker.py"
    code = r'''#!/usr/bin/env python3
"""Self-contained MitoZ annotation worker invoked by Snakemake.

This script handles MitoZ annotation only.  Finalization steps
(ID-rewrite, reorientation, metadata overlay, organism sync) are
executed in the main pipeline process after all workers finish.
"""
import json, re, shlex, shutil, subprocess, sys, time
from pathlib import Path

MITOZ_EXE_NAME = "mitoz"

def which_or_none(exe):
    try: return shutil.which(exe)
    except: return None

def build_tool_cmd(tool_label, exe_name, base_args, conda_env, tool_path, conda_exe="conda"):
    if conda_env:
        cb = which_or_none(conda_exe)
        if not cb: raise FileNotFoundError(f"[{tool_label}] '{conda_exe}' not found")
        return [cb, "run", "-n", conda_env, "--no-capture-output", exe_name] + base_args
    if tool_path:
        p = Path(tool_path).expanduser()
        exe = (p / exe_name).resolve() if p.is_dir() else p.resolve()
        if not exe.exists(): raise FileNotFoundError(f"[{tool_label}] not found: {exe}")
        if exe.suffix.lower() == ".py": return [sys.executable, str(exe)] + base_args
        return [str(exe)] + base_args
    found = which_or_none(exe_name)
    if not found: raise FileNotFoundError(f"[{tool_label}] '{exe_name}' not in PATH")
    return [found] + base_args

def tail_text(s, n=120): return "\n".join((s or "").splitlines()[-n:])

def _is_gb_file(p):
    return p.suffix.lower() in {".gb", ".gbk", ".genbank", ".gbf"}

def find_genbank_like_files(root, preferred_prefix=None):
    out = [p for p in Path(root).rglob("*") if p.is_file() and _is_gb_file(p)]
    if preferred_prefix:
        pref = sorted([p for p in out if p.stem.startswith(preferred_prefix)],
                      key=lambda x: x.stat().st_size if x.exists() else 0, reverse=True)
        rest = sorted([p for p in out if not p.stem.startswith(preferred_prefix)],
                      key=lambda x: x.stat().st_size if x.exists() else 0, reverse=True)
        return pref + rest
    return sorted(out, key=lambda x: x.stat().st_size if x.exists() else 0, reverse=True)

def find_any_summary_txt(d):
    cands = list(Path(d).rglob("summary.txt"))
    if not cands: return None
    # BUG-V0.29-F: keep worker logic consistent with the main helper.  The
    # old len(str(p)) sort could choose a shallow file with a long path.
    # cands.sort(key=lambda p: (len(str(p)), str(p)), reverse=True)
    cands.sort(key=lambda p: (len(p.parts), str(p)), reverse=True)
    return cands[0]

def parse_counts(f):
    txt = Path(f).read_text(encoding="utf-8", errors="ignore")
    m1 = re.search(r"Protein\s+coding\s+genes\s+totally\s+found:\s*(\d+)", txt, re.I)
    m2 = re.search(r"tRNA\s+genes\s+totally\s+found:\s*(\d+)",             txt, re.I)
    m3 = re.search(r"rRNA\s+genes\s+totally\s+found:\s*(\d+)",             txt, re.I)
    # BUG-V0.29-F: old worker returned (0,0,0) on parse miss, preventing the
    # GenBank-count fallback used by the threads backend.
    # return (int(m1.group(1)), int(m2.group(1)), int(m3.group(1))) if m1 and m2 and m3 else (0,0,0)
    return (int(m1.group(1)), int(m2.group(1)), int(m3.group(1))) if m1 and m2 and m3 else None

def counts_from_gb(gb):
    from Bio import SeqIO as BIO
    pcg = trna = rrna = 0
    for rec in BIO.parse(str(gb), "genbank"):
        for feat in rec.features:
            if feat.type == "CDS":    pcg  += 1
            elif feat.type == "tRNA": trna += 1
            elif feat.type == "rRNA": rrna += 1
    return pcg, trna, rrna

def main():
    task = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    args = json.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))
    internal_id   = task["internal_id"]
    phase_dirname = task["phase_dirname"]
    out_root      = Path(task["out_root"]).resolve()
    fasta_file    = Path(task["fasta_file"]).resolve()
    fq1, fq2      = task.get("fq1"), task.get("fq2")
    workdir       = out_root / phase_dirname / internal_id
    workdir.mkdir(parents=True, exist_ok=True)
    logs_dir = workdir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    base_args = [
        "annotate",
        "--workdir",       str(workdir),
        "--outprefix",     internal_id,
        "--thread_number", str(int(args["threads"])),
        "--clade",         str(args["clade"]),
        "--genetic_code",  str(int(args["genetic_code"])),
        "--fastafiles",    str(fasta_file),
    ]
    if args.get("species_name"):
        base_args += ["--species_name", str(args["species_name"])]
    if fq1 and fq2:
        base_args += ["--fq1", str(fq1), "--fq2", str(fq2)]
    for x in args.get("mitoz_extra", []):
        base_args += shlex.split(str(x))

    cmd = build_tool_cmd(
        "MitoZ", MITOZ_EXE_NAME, base_args,
        args.get("conda_env_mitoz"), args.get("mitoz_path"),
        args.get("conda_exe", "conda"))

    t0 = time.time()
    p  = subprocess.run(cmd, cwd=str(workdir),
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    dt = time.time() - t0

    if p.stdout:
        (logs_dir / "mitoz_annotate.stdout.log").write_text(
            p.stdout, encoding="utf-8", errors="ignore")
    if p.stderr:
        (logs_dir / "mitoz_annotate.stderr.log").write_text(
            p.stderr, encoding="utf-8", errors="ignore")
    if p.returncode != 0:
        (logs_dir / "mitoz_annotate.returncode.txt").write_text(
            str(p.returncode) + "\n", encoding="utf-8")
    (logs_dir / "mitoz_annotate.runtime_sec.txt").write_text(
        f"{dt:.2f}\n", encoding="utf-8")

    result = {
        "internal_id": internal_id,
        "sample_name": task.get("sample_name", internal_id),
        "ok": False, "workdir": str(workdir),
        "best_gb": None, "summary_txt": None, "counts": [0,0,0],
        "returncode": int(p.returncode),
    }
    if p.returncode == 0:
        gb_files = find_genbank_like_files(workdir, preferred_prefix=internal_id)
        best_gb  = str(gb_files[0]) if gb_files else None
        summ     = find_any_summary_txt(workdir)
        # BUG-V0.29-F: keep fallback active when summary.txt exists but cannot
        # be parsed.
        # counts   = list(parse_counts(summ)) if summ else None
        parsed   = parse_counts(summ) if summ else None
        counts   = list(parsed) if parsed is not None else None
        if counts is None and best_gb:
            try: counts = list(counts_from_gb(Path(best_gb)))
            except: counts = [0,0,0]
        result.update({"ok": best_gb is not None, "best_gb": best_gb,
                        "summary_txt": str(summ) if summ else None,
                        "counts": counts or [0,0,0]})
    else:
        result["stderr_tail"] = tail_text(p.stderr or "")

    Path(sys.argv[1]).with_suffix(".result.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    Path(sys.argv[1]).with_suffix(".done").write_text("done\n", encoding="utf-8")

if __name__ == "__main__":
    main()
'''
    fp.write_text(code, encoding="utf-8")
    return fp


def render_snakemake_snakefile(
    workdir: Path, task_jsons: List[Path], worker_script: Path, args_json: Path,
    python_exe: str = sys.executable or "python3",
    job_threads: int = 1,
    job_mem_mb: int = 128,
) -> Path:
    snakefile = workdir / "Snakefile"
    rule_lines: List[str] = []
    all_done:   List[str] = []
    for tp in task_jsons:
        done_s  = str(tp.with_suffix(".done"))
        rule_id = re.sub(r"[^A-Za-z0-9_]", "_", tp.stem)
        all_done.append(done_s)
        cmdline = " ".join(shlex.quote(x) for x in [
            python_exe, str(worker_script), str(tp), str(args_json)
        ])
        rule_lines.append(f"""
rule task_{rule_id}:
    input:  task={json.dumps(str(tp), ensure_ascii=False)}
    output: done={json.dumps(done_s, ensure_ascii=False)}
    # BUG-V0.29-G: old code declared threads: 1 while MitoZ receives
    # --thread_number, causing Snakemake to over-schedule CPU.
    # threads: 1
    threads: {max(1, int(job_threads))}
    # BUG-V0.29-G: --resources mem_mb was previously global only; each rule
    # must request memory for Snakemake to enforce the pool.
    resources: mem_mb={int_mem_mb(job_mem_mb)}
    # BUG-V0.29-H: old shell used bare python and unquoted paths:
    # shell:  "python '{{worker_script}}' '{{tp}}' '{{args_json}}'"
    shell:  {json.dumps(cmdline, ensure_ascii=False)}
""")
    header = f"rule all:\n    input: {json.dumps(all_done, ensure_ascii=False)}\n"
    snakefile.write_text(header + "\n".join(rule_lines), encoding="utf-8")
    return snakefile


def run_snakemake_for_annotation_tasks(
    *, tasks: List[AnnotTask], out_root: Path, phase_name: str, args: argparse.Namespace,
) -> Dict[str, Dict]:
    snk_work = out_root / f"_snakemake_{phase_name}"
    snk_work.mkdir(parents=True, exist_ok=True)
    worker_script = render_snakemake_worker_script(snk_work)
    args_json     = snk_work / "args.json"
    args_payload  = {
        "threads": int(args.threads), "clade": args.clade,
        "genetic_code": int(args.genetic_code), "species_name": args.species_name,
        "mitoz_extra": args.mitoz_extra or [], "conda_env_mitoz": args.conda_env_mitoz,
        "mitoz_path": args.mitoz_path, "conda_exe": args.conda_exe,
    }
    args_json.write_text(json.dumps(args_payload, indent=2, ensure_ascii=False) + "\n",
                         encoding="utf-8")
    task_jsons: List[Path] = []
    for t in tasks:
        tp = snk_work / f"{t.internal_id}.json"
        tp.write_text(json.dumps(asdict(t), indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
        task_jsons.append(tp)
    per_job_mem_mb = max(128, int_mem_mb(args.smk_mem_mb) // max(1, int(args.max_tasks)))
    # BUG-V0.29-G/H: pass the real per-job CPU/memory requirements and the
    # current interpreter into the generated Snakefile.
    # snakefile = render_snakemake_snakefile(snk_work, task_jsons, worker_script, args_json)
    snakefile = render_snakemake_snakefile(
        snk_work, task_jsons, worker_script, args_json,
        python_exe=sys.executable or "python3",
        job_threads=int(args.threads),
        job_mem_mb=per_job_mem_mb)
    smk_bin   = args.snakemake_bin or which_or_none("snakemake")
    if not smk_bin:
        raise FileNotFoundError("snakemake executable not found")
    cmd = [
        str(smk_bin), "--snakefile", str(snakefile), "--directory", str(snk_work),
        "--cores", str(max(1, int(args.smk_cores))),
        "--jobs",  str(max(1, int(args.max_tasks))),
        "--resources", f"mem_mb={int_mem_mb(args.smk_mem_mb)}",
        "--rerun-incomplete", "--keep-going", "--printshellcmds",
        "--latency-wait", str(int(args.smk_latency_wait)),
    ]
    smk_extra: List[str] = []
    for token in (getattr(args, "smk_extra", None) or []):
        smk_extra.extend(shlex.split(str(token)))
    if smk_extra:
        logging.info("[SNAKEMAKE_EXTRA] " + " ".join(smk_extra))
        cmd.extend(smk_extra)

    def _run(c: List[str]) -> "subprocess.CompletedProcess[str]":
        logging.info("[SNAKEMAKE] " + " ".join(c))
        return subprocess.run(c, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    p = _run(cmd)
    if p.returncode != 0:
        combined = (p.stdout or "") + (p.stderr or "")
        if "LockException" in combined or "Directory cannot be locked" in combined:
            logging.warning(f"[SNAKEMAKE_LOCK] Stale lock in {snk_work} — attempting --unlock…")
            pu = subprocess.run(
                [str(smk_bin), "--snakefile", str(snakefile),
                 "--directory", str(snk_work), "--unlock"],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if pu.returncode == 0:
                logging.warning("[SNAKEMAKE_LOCK] Unlock OK — retrying…")
                p = _run(cmd)
            else:
                logging.error(f"[SNAKEMAKE_UNLOCK_FAIL] rc={pu.returncode}")

    (snk_work / "snakemake.stdout.log").write_text(p.stdout or "", encoding="utf-8",
                                                    errors="ignore")
    (snk_work / "snakemake.stderr.log").write_text(p.stderr or "", encoding="utf-8",
                                                    errors="ignore")
    if p.returncode != 0:
        raise RuntimeError(
            f"Snakemake failed for phase '{phase_name}'\n"
            f"stdout:\n{tail_text(p.stdout or '')}\n\nstderr:\n{tail_text(p.stderr or '')}")

    results: Dict[str, Dict] = {}
    for tp in task_jsons:
        rp = tp.with_suffix(".result.json")
        if rp.exists():
            obj = json.loads(rp.read_text(encoding="utf-8"))
            results[obj["internal_id"]] = obj
        else:
            iid = tp.stem
            logging.warning(f"[SNAKEMAKE_RESULT_MISSING] {rp} not found for task {iid}")
            results[iid] = {
                "internal_id": iid, "sample_name": iid, "ok": False,
                "workdir": str(out_root / phase_name / iid),
                "best_gb": None, "summary_txt": None, "counts": [0, 0, 0],
            }
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Unified annotation stage dispatcher
# ─────────────────────────────────────────────────────────────────────────────
def run_annotation_stage(
    *, samples: List[SampleInfo], phase_dirname: str, out_root: Path,
    args: argparse.Namespace, fasta_selector,
) -> Dict[str, Dict]:
    infq = Path(args.infq).resolve() if getattr(args, "infq", None) else None
    suffix1, suffix2 = parse_suffix_fq(args.suffix_fq)
    backend = getattr(args, "scheduler_backend", "threads")
    if backend == "auto":
        backend = "snakemake" if snakemake_available(getattr(args, "snakemake_bin", None)) \
                  else "threads"
    logging.info(f"[ANNOT_STAGE] phase={phase_dirname}  backend={backend}  n={len(samples)}")

    if backend == "snakemake":
        tasks: List[AnnotTask] = []
        for si in samples:
            fq1 = fq2 = None
            if infq:
                fq1p, fq2p = find_fq_pair_for_sample(
                    si.sample_name, infq, args.fq_position, suffix1, suffix2)
                fq1 = str(fq1p) if fq1p else None
                fq2 = str(fq2p) if fq2p else None
            tasks.append(AnnotTask(
                internal_id=si.internal_id, sample_name=si.sample_name,
                fasta_file=str(fasta_selector(si)), phase_dirname=phase_dirname,
                out_root=str(out_root), fq1=fq1, fq2=fq2))
        return run_snakemake_for_annotation_tasks(
            tasks=tasks, out_root=out_root, phase_name=phase_dirname, args=args)

    results: Dict[str, Dict] = {}
    total = len(samples)
    done_counter: List = [0, threading.Lock()]

    def _task(si: SampleInfo) -> Tuple:
        fq1 = fq2 = None
        if infq:
            fq1, fq2 = find_fq_pair_for_sample(
                si.sample_name, infq, args.fq_position, suffix1, suffix2)
        return mitoz_annotate_one(
            si, fasta_selector(si), out_root, phase_dirname, args, fq1, fq2,
            total=total, done_counter=done_counter)

    with cf.ThreadPoolExecutor(max_workers=int(args.max_tasks)) as ex:
        futs = {ex.submit(_task, si): si for si in samples}
        for fut in cf.as_completed(futs):
            si = futs[fut]
            try:
                internal_id, ok, workdir, best_gb, summary_txt, counts = fut.result()
            except Exception as exc:
                logging.error(f"[TASK_EXCEPTION] {si.sample_name}: {exc}")
                internal_id, ok, workdir, best_gb, summary_txt, counts = (
                    si.internal_id, False,
                    out_root / phase_dirname / si.internal_id,
                    None, None, (0, 0, 0))
            results[internal_id] = {
                "internal_id": internal_id, "sample_name": si.sample_name,
                "ok": ok, "workdir": str(workdir),
                "best_gb":      str(best_gb)     if best_gb     else None,
                "summary_txt":  str(summary_txt) if summary_txt else None,
                "counts": counts,
            }
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Final GenBank output
# ─────────────────────────────────────────────────────────────────────────────
def finalize_genbank_output(
    *,
    sample_name:       str,
    body_src_gb:       Path,
    metadata_src_gb:   Optional[Path],
    dest_gb:           Path,
    which_gene_first:  str,
    force_forward:     bool,
    preserve_metadata: bool,
    smart_metadata:    bool = False,
    auto_accept:       bool = False,
) -> Tuple[int, int, int, str, bool]:
    """Produce the final, publication-grade GenBank file for one sample.

    Steps
    -----
    1. ID rewrite      — replace internal S0xxxx IDs with sample_name.
    2. Reorientation   — guarantee the record starts at which_gene_first.
    3. Metadata overlay — (when preserve_metadata) replace GB header section.
    3.5 Organism sync  — fix /organism= inside the FEATURES source feature.
    4. Recount         — recount PCG / tRNA / rRNA from the definitive file.

    BUG-NEW-F (v0.27): bare shutil.copy2 calls replaced with _safe_copy.
    """
    pcg = trna = rrna = 0
    orientation_ok = False

    td = Path(tempfile.mkdtemp(prefix="mitoz_final_"))
    dest_gb.parent.mkdir(parents=True, exist_ok=True)
    draft_gb = dest_gb.with_name(f".{dest_gb.name}.tmp")
    committed = False
    try:
        if draft_gb.exists():
            draft_gb.unlink()

        # ── Step 1: ID rewrite ────────────────────────────────────────
        id_gb = td / "id_rewrite.gbf"
        try:
            pcg, trna, rrna = rewrite_genbank_record_ids(body_src_gb, id_gb, sample_name)
        except Exception as exc:
            raise RuntimeError(f"ID rewrite failed for {sample_name}: {exc}") from exc

        # ── Step 2: Final guaranteed reorientation ────────────────────
        reor_gb    = td / "reoriented.gbf"
        no_gene_gb = td / "no_gene.gbf"
        try:
            n_ok, n_no, _ = reorient_one_gb(
                id_gb, reor_gb, no_gene_gb, which_gene_first, force_forward)
            if n_ok > 0 and reor_gb.exists():
                orientation_ok = True
                if n_no > 0 and no_gene_gb.exists():
                    merged_reor = td / "merged_reoriented.gbf"
                    merge_gb_files([reor_gb, no_gene_gb], merged_reor)
                    working_gb = merged_reor
                    logging.info(
                        f"[FINAL_MULTICONTIG] {sample_name}: merged {n_ok} reoriented "
                        f"+ {n_no} un-reoriented contig(s) into final GB.")
                else:
                    working_gb = reor_gb
            else:
                working_gb = id_gb
                logging.info(
                    f"[FINAL_NO_GENE] {sample_name}: '{which_gene_first}' not found "
                    f"in chosen_gb; orientation unchanged.")
        except Exception as exc:
            logging.warning(f"[FINAL_REORIENT_FAIL] {sample_name}: {exc} — using ID-rewritten GB")
            working_gb = id_gb

        # ── Step 3: Optional metadata overlay (header section, text-level) ──
        # BUG-V0.29-L: old code copied each intermediate directly to dest_gb.
        # If a later step failed, a half-final file could remain at the final
        # path.  Write a same-directory draft and replace dest_gb at the end.
        # dest_gb.parent.mkdir(parents=True, exist_ok=True)
        metadata_mode = "rewritten_no_metadata"

        if preserve_metadata and metadata_src_gb and metadata_src_gb.exists():
            if gb_file_has_meaningful_metadata(metadata_src_gb):
                meta_out = td / "with_metadata.gbf"
                ok, n_replaced, n_unmatched, unmatched, reason, compat = transfer_gb_metadata_file(
                    metadata_src_gb=metadata_src_gb,
                    new_body_gb=working_gb,
                    output_gb=meta_out,
                    smart=smart_metadata,
                    auto_accept=auto_accept,
                )
                if ok and n_replaced > 0:
                    # _safe_copy(meta_out, dest_gb)          # BUG-NEW-F
                    _safe_copy(meta_out, draft_gb)         # BUG-V0.29-L
                    metadata_mode = "metadata_transferred"
                    if compat:
                        metadata_mode += "+compat"
                    if unmatched:
                        logging.warning(
                            f"[META_UNMATCHED] {sample_name}: {len(unmatched)} record(s) "
                            f"unmatched: {unmatched}")
                elif ok and n_replaced == 0:
                    # _safe_copy(working_gb, dest_gb)        # BUG-NEW-F
                    _safe_copy(working_gb, draft_gb)       # BUG-V0.29-L
                    metadata_mode = "metadata_no_match"
                    logging.info(
                        f"[META_NO_MATCH] {sample_name}: 0 records matched in "
                        f"{metadata_src_gb.name}; using ID-rewritten GB.")
                else:
                    # _safe_copy(working_gb, dest_gb)        # BUG-NEW-F
                    _safe_copy(working_gb, draft_gb)       # BUG-V0.29-L
                    metadata_mode = "metadata_failed"
                    logging.warning(f"[META_FAIL] {sample_name}: {reason}")
            else:
                # _safe_copy(working_gb, dest_gb)            # BUG-NEW-F
                _safe_copy(working_gb, draft_gb)           # BUG-V0.29-L
                metadata_mode = "metadata_src_not_usable"
        else:
            # _safe_copy(working_gb, dest_gb)                # BUG-NEW-F
            _safe_copy(working_gb, draft_gb)               # BUG-V0.29-L

        # ── Step 3.5: Sync /organism= in FEATURES source (BUG-12 fix) ──
        try:
            org_fixed = td / "org_fixed.gbf"
            applied = False

            if metadata_src_gb and metadata_src_gb.exists():
                org_index = build_organism_index([metadata_src_gb])
                if org_index:
                    n_rep, n_nf, n_tot = replace_organism_in_gb_file(
                        # dest_gb, org_fixed, org_index)
                        draft_gb, org_fixed, org_index)
                    if n_rep > 0:
                        # _safe_copy(org_fixed, dest_gb)     # BUG-NEW-F
                        _safe_copy(org_fixed, draft_gb)    # BUG-V0.29-L
                        applied = True
                        logging.info(
                            f"[ORG_SYNC] {sample_name}: /organism= restored in "
                            f"{n_rep}/{n_tot} record(s) from original GB"
                            + (f" ({n_nf} not matched)" if n_nf else ""))

            if not applied:
                # n_synced = _sync_organism_qualifier(dest_gb, org_fixed)
                n_synced = _sync_organism_qualifier(draft_gb, org_fixed)
                if n_synced > 0:
                    # _safe_copy(org_fixed, dest_gb)         # BUG-NEW-F
                    _safe_copy(org_fixed, draft_gb)        # BUG-V0.29-L
                    logging.info(
                        f"[ORG_SYNC] {sample_name}: /organism= synced from header "
                        f"in {n_synced} record(s)")

        except Exception as exc:
            logging.warning(
                f"[ORG_SYNC_FAIL] {sample_name}: organism qualifier sync failed: {exc}")

        # ── Step 4: recount from definitive file ──────────────────────
        try:
            # pcg, trna, rrna = counts_from_genbank(dest_gb)
            pcg, trna, rrna = counts_from_genbank(draft_gb)
        except Exception:
            pass
        os.replace(str(draft_gb), str(dest_gb))
        committed = True

    finally:
        if not committed:
            try:
                draft_gb.unlink()
            except Exception:
                pass
        shutil.rmtree(td, ignore_errors=True)

    return pcg, trna, rrna, metadata_mode, orientation_ok


# ─────────────────────────────────────────────────────────────────────────────
# Workflows
# ─────────────────────────────────────────────────────────────────────────────
def workflow_annotate(args: argparse.Namespace) -> int:
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "pipeline.log")
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))
    preflight(args)

    in_path     = Path(args.input).resolve()
    split_dir   = out_root / "00.internal_fasta"
    gb_conv_dir = out_root / "00.gb_converted_fasta"

    in_files, input_type, acc_to_orig_gb = iter_inputs(in_path, gb_to_fasta_dir=gb_conv_dir)
    logging.info(f"[INPUT_TYPE] {input_type}  ({len(in_files)} file(s))")
    samples, id_map = split_or_sanitize_fastas(
        in_files, split_dir, args.if_split_fasta == "yes", int(args.id_maxlen),
        acc_to_orig_gb=acc_to_orig_gb)
    logging.info(f"[ID_MAP] {id_map}  (total: {len(samples)})")

    results = run_annotation_stage(
        samples=samples, phase_dirname="01.mitoz_anno1",
        out_root=out_root, args=args,
        fasta_selector=lambda si: si.internal_fasta)

    summary    = RunSummary()
    final_dir  = out_root / "04.final_gb"
    final_dir.mkdir(parents=True, exist_ok=True)
    summary_rows: List[str] = []

    for si in samples:
        r       = results.get(si.internal_id, {})
        best_gb = Path(r["best_gb"]) if r.get("best_gb") else None
        if not r.get("ok") or not best_gb or not best_gb.exists():
            summary.annotation_failed.append(si.sample_name); continue
        dest = final_dir / f"{si.sample_name}.gbf"
        try:
            pcg, trna, rrna, mode, orient_ok = finalize_genbank_output(
                sample_name=si.sample_name, body_src_gb=best_gb,
                metadata_src_gb=si.orig_gb_path, dest_gb=dest,
                which_gene_first=args.which_gene_first,
                force_forward=args.force_forward,
                preserve_metadata=(args.preserve_phase1_metadata == "yes")
                    if hasattr(args, "preserve_phase1_metadata") else False,
            )
            summary_rows.append(
                f"{si.sample_name}\t{pcg}\t{trna}\t{rrna}\t{'yes' if orient_ok else 'no'}\t{mode}")
        except Exception as exc:
            logging.error(f"[FINAL_WRITE_FAIL] {si.sample_name}: {exc}")
            summary.final_write_failed.append(si.sample_name)

    write_nonempty_list(out_root / "annotation_failed.txt",
                        summary.annotation_failed, "# MitoZ non-zero exit")
    write_nonempty_list(out_root / "final_write_failed.txt",
                        summary.final_write_failed, "# Failed during final write")
    write_nonempty_table(
        out_root / "annotation-summary.tsv",
        "sample_name\tPCG\ttRNA\trRNA\tcox1_orientation_ok\tmetadata_mode",
        summary_rows)
    # BUG-NEW-D (v0.27): use rglob consistently with workflow_run so that
    # no_{gene}_gene/ subdirectory records are also included in the merged
    # archive.  Use final_dir.glob("*.gbf") if you want only the top-level
    # (successfully-oriented) records.
    all_final_gbs = sorted(
        p for p in final_dir.rglob("*.gbf")
        if p.is_file() and p.stat().st_size > 0)
    if all_final_gbs:
        merge_gb_files(all_final_gbs, out_root / "all.final.gbf")
    (out_root / "run_params.json").write_text(
        json.dumps(vars(args), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    summary.log()
    logging.info("[DONE]")
    return 0


def workflow_reorient(args: argparse.Namespace) -> int:
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "pipeline.log")
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))

    in_path  = Path(args.input).resolve()
    gb_files = collect_genbank_files(in_path)
    if not gb_files:
        logging.error("[ERROR] No GB files found.")
        return 2

    out_dir     = out_root / "reoriented"
    no_gene_dir = out_root / f"no_{args.which_gene_first}_gene"
    out_dir.mkdir(parents=True, exist_ok=True)
    no_gene_dir.mkdir(parents=True, exist_ok=True)

    summary      = RunSummary()
    summary_rows: List[str] = []

    for gb_in in gb_files:
        sample     = safe_name(gb_in.stem)
        gb_out     = out_dir     / f"{sample}.gbf"
        gb_no_gene = no_gene_dir / f"{sample}.gbf"
        try:
            n_ok, n_no, no_gene_ids = reorient_one_gb(
                gb_in, gb_out, gb_no_gene, args.which_gene_first, args.force_forward)
            if n_no > 0:
                logging.warning(f"[NO_TARGET_GENE] '{sample}': "
                                 f"{n_no} record(s) missing '{args.which_gene_first}'")
                summary.no_target_gene.append(sample)
            if n_ok > 0 and gb_out.exists():
                pcg, trna, rrna = counts_from_genbank(gb_out)
                summary_rows.append(f"{sample}\t{pcg}\t{trna}\t{rrna}\tyes")
            elif n_ok == 0 and n_no > 0:
                _safe_copy(gb_in, gb_out)
                pcg, trna, rrna = counts_from_genbank(gb_out)
                summary_rows.append(f"{sample}\t{pcg}\t{trna}\t{rrna}\tno")
                logging.warning(f"[FALLBACK_COPY] '{sample}': original copied unchanged")
        except Exception as exc:
            logging.error(f"[FAIL] {sample}: {exc}")
            # BUG-V0.28-F: use reorient_failed — annotation_failed is
            # semantically incorrect here (no MitoZ annotation is performed
            # in the reorient workflow).
            summary.reorient_failed.append(sample)

    write_nonempty_list(out_root / "no_target_gene_samples.txt",
                        summary.no_target_gene, f"# samples without {args.which_gene_first}")
    write_nonempty_list(out_root / "failed_samples.txt",
                        summary.all_problem_samples(), "# all problematic samples")
    write_nonempty_table(out_root / "annotation-summary.tsv",
                         "sample_name\tPCG\ttRNA\trRNA\tcox1_orientation_ok",
                         summary_rows)
    summary.log(); logging.info("[DONE]")
    return 0


def workflow_run(args: argparse.Namespace) -> int:
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "pipeline.log")
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))
    preflight(args)

    in_path     = Path(args.input).resolve()
    split_dir   = out_root / "00.internal_fasta"
    gb_conv_dir = out_root / "00.gb_converted_fasta"

    in_files, input_type, acc_to_orig_gb = iter_inputs(in_path, gb_to_fasta_dir=gb_conv_dir)
    logging.info(f"[INPUT_TYPE] {input_type}  ({len(in_files)} file(s))")
    samples, id_map = split_or_sanitize_fastas(
        in_files, split_dir, args.if_split_fasta == "yes", int(args.id_maxlen),
        acc_to_orig_gb=acc_to_orig_gb)
    total = len(samples)
    logging.info(f"[ID_MAP] {id_map}  (total: {total})")

    summary = RunSummary()

    results1 = run_annotation_stage(
        samples=samples, phase_dirname="01.mitoz_anno1",
        out_root=out_root, args=args,
        fasta_selector=lambda si: si.internal_fasta)
    for si in samples:
        if not results1.get(si.internal_id, {}).get("ok"):
            summary.annotation_failed.append(si.sample_name)

    logging.info(
        f"[PHASE1b] reorient gene='{args.which_gene_first}' force_forward={args.force_forward}")
    reor_fasta_dir  = out_root / "00.reoriented_fasta_internal"
    reor_fasta_dir.mkdir(parents=True, exist_ok=True)
    reor_ok_ids:    List[str]        = []
    no_gene_iids:   set              = set()
    # BUG-V0.29-C: keep partial no-gene contigs so they can be merged back
    # into the final body.  The old workflow used only gb_reor downstream and
    # silently dropped these records.
    partial_no_gene_gb: Dict[str, Path] = {}

    for si in samples:
        r1      = results1.get(si.internal_id, {})
        best_gb = Path(r1["best_gb"]) if r1.get("best_gb") else None
        if not best_gb or not best_gb.exists():
            continue
        reor_dir   = out_root / "01.mitoz_anno1" / si.internal_id / "reoriented"
        reor_dir.mkdir(parents=True, exist_ok=True)
        gb_reor    = reor_dir / f"{si.internal_id}.reoriented.gbf"
        gb_no_gene = reor_dir / f"{si.internal_id}.no_{args.which_gene_first}.gbf"
        try:
            n_ok, n_no, no_gene_ids = reorient_one_gb(
                best_gb, gb_reor, gb_no_gene, args.which_gene_first, args.force_forward)
            if n_ok == 0 and n_no > 0:
                logging.warning(
                    f"[NO_TARGET_GENE] '{si.sample_name}' ({si.internal_id}): "
                    f"ALL {n_no} record(s) missing '{args.which_gene_first}' — skipping Phase-2")
                summary.no_target_gene.append(si.sample_name)
                no_gene_iids.add(si.internal_id)
            elif n_no > 0:
                logging.info(
                    f"[PARTIAL_NO_GENE] '{si.sample_name}' ({si.internal_id}): "
                    f"{n_no} minor contig(s) lack '{args.which_gene_first}'; "
                    f"primary reoriented OK — proceeding to Phase-2.")
                if gb_no_gene.exists():
                    partial_no_gene_gb[si.internal_id] = gb_no_gene
            if n_ok > 0 and gb_reor.exists():
                reor_ok_ids.append(si.internal_id)
                gb_to_internal_fasta(gb_reor, reor_fasta_dir / f"{si.internal_id}.fasta",
                                     si.internal_id)
        except Exception as exc:
            logging.error(f"[REORIENT_FAIL] {si.sample_name} ({si.internal_id}): {exc}")
            # BUG-V0.29-D: a reorientation exception is not the same as a
            # biological no-target-gene case.  Keep the old line commented for
            # traceability and record it in reorient_failed instead.
            # summary.no_target_gene.append(si.sample_name)
            summary.reorient_failed.append(si.sample_name)
            no_gene_iids.add(si.internal_id)

    do_phase2 = (args.if_reannotation_2 == "yes")
    results2:  Dict[str, Dict] = {}
    if do_phase2:
        phase2_samples = [si for si in samples if si.internal_id in reor_ok_ids]
        logging.info(f"[PHASE2] re-annotation: {len(phase2_samples)} samples")
        results2 = run_annotation_stage(
            samples=phase2_samples, phase_dirname="03.mitoz_anno2_reoriented",
            out_root=out_root, args=args,
            fasta_selector=lambda si: reor_fasta_dir / f"{si.internal_id}.fasta")
        for si in phase2_samples:
            if not results2.get(si.internal_id, {}).get("ok"):
                logging.warning(
                    f"[PHASE2_FAIL_FALLBACK] '{si.sample_name}': "
                    f"phase2 failed, fallback to phase1b reoriented GB")

    final_dir         = out_root / "04.final_gb"
    no_gene_final_dir = final_dir / f"no_{args.which_gene_first}_gene"
    final_dir.mkdir(parents=True, exist_ok=True)
    no_gene_final_dir.mkdir(parents=True, exist_ok=True)

    summary_rows:    List[str] = []
    meta_log_rows:   List[str] = []

    for si in samples:
        sample      = si.sample_name
        iid         = si.internal_id
        is_no_gene  = iid in no_gene_iids
        r1          = results1.get(iid, {})
        metadata_src = si.orig_gb_path

        chosen_gb: Optional[Path] = None
        src_tag   = "phase1_raw"

        if do_phase2 and not is_no_gene:
            r2 = results2.get(iid, {})
            if r2.get("best_gb") and Path(r2["best_gb"]).exists():
                chosen_gb = Path(r2["best_gb"])
                src_tag   = "phase2_mitoz"

        if chosen_gb is None:
            gb_reor = (out_root / "01.mitoz_anno1" / iid / "reoriented"
                       / f"{iid}.reoriented.gbf")
            if gb_reor.exists():
                chosen_gb = gb_reor
                src_tag   = "phase1b_reoriented"

        if chosen_gb is None and r1.get("best_gb") and Path(r1["best_gb"]).exists():
            chosen_gb = Path(r1["best_gb"])

        if chosen_gb is None or not chosen_gb.exists():
            if sample not in summary.annotation_failed:
                summary.final_write_failed.append(sample)
            continue

        body_for_final = chosen_gb
        extra_no_gene_gb = partial_no_gene_gb.get(iid)
        if extra_no_gene_gb and extra_no_gene_gb.exists():
            chosen_is_phase1_raw = False
            if r1.get("best_gb"):
                try:
                    chosen_is_phase1_raw = chosen_gb.resolve() == Path(r1["best_gb"]).resolve()
                except Exception:
                    chosen_is_phase1_raw = False
            if not chosen_is_phase1_raw:
                merged_body_dir = out_root / "04.final_body_merged"
                merged_body_dir.mkdir(parents=True, exist_ok=True)
                merged_body = merged_body_dir / f"{iid}.with_no_gene.gbf"
                merge_gb_files([chosen_gb, extra_no_gene_gb], merged_body)
                body_for_final = merged_body
                src_tag += "+partial_no_gene_merged"
                logging.info(
                    f"[FINAL_PARTIAL_NO_GENE_MERGE] '{sample}': merged "
                    f"{extra_no_gene_gb.name} back into final body.")

        dest = (no_gene_final_dir if is_no_gene else final_dir) / f"{sample}.gbf"

        try:
            pcg, trna, rrna, meta_mode, orient_ok = finalize_genbank_output(
                sample_name=sample,
                # BUG-V0.29-C: old code passed chosen_gb directly, dropping
                # partial no-gene contigs created during Phase-1b.
                # body_src_gb=chosen_gb,
                body_src_gb=body_for_final,
                metadata_src_gb=metadata_src,
                dest_gb=dest,
                which_gene_first=args.which_gene_first,
                force_forward=args.force_forward,
                preserve_metadata=(args.preserve_phase1_metadata == "yes"),
                smart_metadata=(args.smart_metadata == "yes"),
                auto_accept=(args.auto_accept_metadata == "yes"),
            )
            summary_rows.append(
                f"{sample}\t{pcg}\t{trna}\t{rrna}\t"
                f"{'yes' if orient_ok else 'no'}\t{src_tag}\t{meta_mode}")
            meta_log_rows.append(
                # BUG-V0.29-C: log the actual final body source after partial
                # no-gene contigs are merged back.
                # f"{sample}\t{iid}\t{metadata_src or 'NA'}\t{chosen_gb}\t{dest}\t{meta_mode}")
                f"{sample}\t{iid}\t{metadata_src or 'NA'}\t{body_for_final}\t{dest}\t{meta_mode}")

            if "transferred" in meta_mode:
                summary.metadata_transferred.append(sample)
                if "compat" in meta_mode:
                    summary.metadata_compat_match.append(sample)
                    # BUG-V0.29-I: old code counted compat matches but never
                    # filled metadata_key_mismatch, so the audit file stayed
                    # empty even when metadata keys differed.
                    summary.metadata_key_mismatch.append(sample)
            elif meta_mode in ("metadata_no_match", "metadata_src_not_usable"):
                summary.metadata_skipped.append(sample)
            elif meta_mode == "metadata_failed":
                summary.metadata_transfer_failed.append(sample)

        except Exception as exc:
            logging.error(f"[FINAL_WRITE_FAIL] {sample}: {exc}")
            summary.final_write_failed.append(sample)

    write_nonempty_table(
        out_root / "annotation-summary.tsv",
        "sample_name\tPCG\ttRNA\trRNA\tcox1_orientation_ok\tbody_source\tmetadata_mode",
        summary_rows)
    write_nonempty_table(
        out_root / "metadata_transfer_log.tsv",
        "sample_name\tinternal_id\tmetadata_src\tbody_src\tdest\tmode",
        meta_log_rows)
    write_nonempty_list(out_root / "annotation_failed.txt",
                        summary.annotation_failed, "# MitoZ annotation returned non-zero exit")
    write_nonempty_list(out_root / f"no_{args.which_gene_first}_samples.txt",
                        summary.no_target_gene,
                        f"# samples where '{args.which_gene_first}' was not found")
    write_nonempty_list(out_root / "final_write_failed.txt",
                        summary.final_write_failed, "# failed during final GB write")
    write_nonempty_list(out_root / "reorient_failed.txt",
                        summary.reorient_failed, "# failed during reorientation")
    write_nonempty_list(out_root / "metadata_transfer_failed.txt",
                        summary.metadata_transfer_failed, "# metadata transfer failed")
    write_nonempty_list(out_root / "metadata_key_mismatch.txt",
                        summary.metadata_key_mismatch,
                        "# accession key mismatch during metadata transfer")
    write_nonempty_list(out_root / "all_problem_samples.txt",
                        summary.all_problem_samples(), "# all problematic samples")
    (out_root / "run_params.json").write_text(
        json.dumps(vars(args), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    # Merge all final GBs including no_{gene}_gene/ subdirectory.
    # Use final_dir.glob("*.gbf") if you want only successfully-oriented records.
    all_final_gbs = sorted(
        p for p in final_dir.rglob("*.gbf")
        if p.is_file() and p.stat().st_size > 0
    )
    if all_final_gbs:
        n_merged = merge_gb_files(all_final_gbs, out_root / "all.final.gbf")
        logging.info(f"[MERGED_ALL] {out_root / 'all.final.gbf'}  ({n_merged} records)")
    else:
        logging.warning("[MERGED_ALL] No final GB files to merge.")

    n_ok = len(summary_rows)
    logging.info(f"[FINAL_DIR] {final_dir}")
    logging.info(f"[OUTPUT_COUNT] {n_ok}/{total}")
    logging.info(f"[METADATA_MODE] preserve={args.preserve_phase1_metadata}")
    summary.log()
    logging.info("[DONE]")
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# Standalone transfer_metadata workflow
# ─────────────────────────────────────────────────────────────────────────────
def workflow_transfer_metadata(args: argparse.Namespace) -> int:
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "transfer_metadata.log")
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))

    metadata_files = collect_genbank_files(Path(args.metadata).resolve())
    body_files     = collect_genbank_files(Path(args.new_body).resolve())
    if not metadata_files:
        logging.error(f"[ERROR] No metadata GB files: {args.metadata}"); return 2
    if not body_files:
        logging.error(f"[ERROR] No body GB files: {args.new_body}"); return 2

    out_dir = out_root / "transferred_gb"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_rows: List[str] = []

    def _process(body_file: Path, meta_file: Optional[Path], out_file: Path) -> str:
        if args.if_replace_metadata == "no":
            if args.overwrite_new_body == "no":
                _safe_copy(body_file, out_file)
            log_rows.append(f"{body_file}\t{meta_file or 'NA'}\t{out_file}\tskipped_by_user\t0\t0")
            return "skip"
        if meta_file is None:
            if args.overwrite_new_body == "no":
                _safe_copy(body_file, out_file)
            log_rows.append(f"{body_file}\tNA\t{out_file}\tno_match\t0\t0")
            return "no_match"
        ok, n_rep, n_un, unmatched, reason, compat = transfer_gb_metadata_file(
            meta_file, body_file, out_file,
            smart=(args.smart_metadata == "yes"), auto_accept=(args.auto_accept_metadata == "yes"))
        if not ok and args.overwrite_new_body == "no":
            _safe_copy(body_file, out_file)
        log_rows.append(
            f"{body_file}\t{meta_file}\t{out_file}\t"
            f"{'ok' if ok else 'no_usable_metadata'}\t{n_rep}\t{n_un}")
        return "ok" if ok else "skip"

    is_single = (len(metadata_files) == 1 and len(body_files) == 1
                 and Path(args.metadata).is_file() and Path(args.new_body).is_file())
    if is_single:
        out_name = getattr(args, "output_name", None) or body_files[0].name
        out_file = body_files[0] if args.overwrite_new_body == "yes" else (out_dir / out_name)
        _process(body_files[0], metadata_files[0], out_file)
        write_nonempty_table(out_root / "transfer_metadata.tsv",
                             "body_file\tmetadata_file\toutput_file\tstatus\treplaced\tunmatched",
                             log_rows)
        return 0

    # BUG-V0.29-J: old dict comprehension silently overwrote duplicate stems
    # from different subdirectories, causing body files to receive the wrong
    # header.  Keep the old line commented for auditability.
    # meta_by_stem = {p.stem: p for p in metadata_files}
    meta_by_stem: Dict[str, Path] = {}
    dup_stems: Dict[str, List[Path]] = {}
    for mf in metadata_files:
        if mf.stem in meta_by_stem:
            dup_stems.setdefault(mf.stem, [meta_by_stem[mf.stem]]).append(mf)
        else:
            meta_by_stem[mf.stem] = mf
    if dup_stems:
        for stem, paths in sorted(dup_stems.items()):
            logging.error(
                f"[DUP_METADATA_STEM] '{stem}' appears in multiple metadata files: "
                + ", ".join(str(p) for p in paths))
        return 2
    n_ok = n_skip = n_fail = 0
    for bf in body_files:
        matched  = meta_by_stem.get(bf.stem)
        out_file = bf if args.overwrite_new_body == "yes" else (out_dir / bf.name)
        try:
            status = _process(bf, matched, out_file)
            if status == "ok": n_ok += 1
            else:              n_skip += 1
        except Exception as exc:
            if args.overwrite_new_body == "no":
                _safe_copy(bf, out_file)
            n_fail += 1
            log_rows.append(f"{bf}\t{matched or 'NA'}\t{out_file}\tfail\t0\t0")
    write_nonempty_table(out_root / "transfer_metadata.tsv",
                         "body_file\tmetadata_file\toutput_file\tstatus\treplaced\tunmatched",
                         log_rows)
    logging.info(f"[SUMMARY] transferred={n_ok}  skipped={n_skip}  failed={n_fail}")
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# Standalone replace_organism workflow  (v0.25 / v0.26 / v0.27)
# ─────────────────────────────────────────────────────────────────────────────
def workflow_replace_organism(args: argparse.Namespace) -> int:
    """Replace /organism= qualifiers in FEATURES source blocks.

    See module docstring and argument help for full description.
    When using the 'run' subcommand this fix is applied automatically as
    Step 3.5 of finalize_genbank_output.
    """
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "replace_organism.log")
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))

    meta_files = collect_genbank_files(Path(args.metadata).resolve())
    body_files = collect_genbank_files(Path(args.new_body).resolve())
    if not meta_files:
        logging.error(f"[ERROR] No metadata GB files: {args.metadata}"); return 2
    if not body_files:
        logging.error(f"[ERROR] No body GB files: {args.new_body}"); return 2

    logging.info(f"[ORG_INDEX] Building organism+source index from {len(meta_files)} file(s)…")
    org_index = build_organism_index(meta_files)
    logging.info(f"[ORG_INDEX] {len(org_index)} accession variant(s) indexed.")
    if not org_index:
        logging.error("[ERROR] Organism index is empty — check metadata GB files."); return 2

    out_dir = out_root / "replaced_organism"
    out_dir.mkdir(parents=True, exist_ok=True)

    log_rows: List[str] = []
    n_ok = n_miss = n_fail = 0

    for bf in body_files:
        out_file = bf if args.overwrite_new_body == "yes" else (out_dir / bf.name)
        try:
            n_rep, n_nf, n_tot = replace_organism_in_gb_file(bf, out_file, org_index)
            # BUG-V0.28-B: append empty note column so all rows have 6 fields
            log_rows.append(f"{bf.name}\t{out_file}\t{n_rep}\t{n_nf}\t{n_tot}\t")
            if n_rep > 0:
                n_ok += 1
                logging.info(
                    f"[DONE] {bf.name}: {n_rep}/{n_tot} record(s) updated, "
                    f"{n_nf} unchanged.")
            else:
                n_miss += 1
                logging.warning(
                    f"[NO_MATCH] {bf.name}: 0 records matched — organism/source not replaced.")
        except Exception as exc:
            n_fail += 1
            logging.error(f"[FAIL] {bf.name}: {exc}")
            # BUG-V0.28-B: FAIL message goes in the dedicated note column (col 6),
            # not appended to the total column (col 5) which must remain numeric.
            log_rows.append(f"{bf.name}\t{out_file}\t0\t0\t0\tFAIL: {exc}")
            _safe_copy(bf, out_file)

    write_nonempty_table(
        out_root / "replace_organism.tsv",
        # BUG-V0.28-B: header now declares 6 columns, matching all data rows
        "source_file\toutput_file\treplaced\tnot_found\ttotal\tnote",
        log_rows)
    logging.info(
        f"[SUMMARY] files_updated={n_ok}  no_match={n_miss}  failed={n_fail}")
    logging.info("[DONE]")
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# Argument parser
# ─────────────────────────────────────────────────────────────────────────────
_CLADE_DEFAULT = "Annelida-segmented-worms"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="batch_mitoz.py",
        description=(
            "Integrated MitoZ batch annotation pipeline — publication-grade\n\n"
            "  annotate          — batch MitoZ annotation\n"
            "  reorient          — reorient .gbf to start at a target gene\n"
            "  run               — full pipeline (annotate → reorient → reanno → final GB)\n"
            "  transfer_metadata — standalone header-transfer between .gb files\n"
            "  replace_organism  — restore /organism= in FEATURES source from original GB\n\n"
            "Snakemake support\n"
            "  Use --scheduler_backend snakemake (or auto, which auto-detects).\n"
            "  Snakemake parallelises the MitoZ annotation phase; all post-\n"
            "  processing (reorientation, metadata overlay, organism sync) runs\n"
            "  in the main process, so all bug-fixes apply regardless of backend.\n"
            "  Cluster submission: --smk_extra '--profile slurm'\n\n"
            "Output files (non-empty only):\n"
            "  04.final_gb/*.gbf         — guaranteed cox1-first records\n"
            "  annotation-summary.tsv    — PCG / tRNA / rRNA counts (final GB)\n"
            "  metadata_transfer_log.tsv — metadata substitution audit log\n"
            "  all.final.gbf              — merged archive of all final records\n"
            "  *_failed.txt / *_mismatch.txt — problem-sample lists\n\n"
            "Metadata preservation (v0.27)\n"
            "  When GenBank files are supplied as input and\n"
            "  --preserve_phase1_metadata yes (default for 'run'), the final\n"
            "  output records carry the original header (SOURCE, ORGANISM,\n"
            "  DEFINITION, REFERENCE …) AND have the /organism= qualifier\n"
            "  inside FEATURES source restored from the original file.\n"
            "  MitoZ placeholder strings (e.g. 'Test sp.') are eliminated.\n\n"
            "Citation:\n"
            "  MitoZ:      Meng et al. (2019) NAR doi:10.1093/nar/gkz173\n"
            "  Snakemake:  Mölder et al. (2021) F1000Res doi:10.12688/f1000research.29032.2\n"
            "  Biopython:  Cock et al. (2009) Bioinformatics doi:10.1093/bioinformatics/btp163\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    def add_common(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("-i", "--input",  required=True,
                        help="Input FASTA file/dir or GenBank file/dir.")
        sp.add_argument("-o", "--out_root", default="mitoz_batch_out",
                        help="Output root directory (default: mitoz_batch_out).")
        sp.add_argument("--if_split_fasta", choices=["yes", "no"], default="yes")
        # BUG-V0.29-E: old parser accepted id_maxlen < 6, which makes S00001
        # style internal IDs collide after truncation.
        # sp.add_argument("--id_maxlen", type=int, default=15,
        sp.add_argument("--id_maxlen", type=id_maxlen_arg, default=15,
                        help="Max length of internal sample IDs (default: 15).")
        sp.add_argument("--which_gene_first", default="cox1",
                        help="Gene used as sequence start (default: cox1).")
        sp.add_argument("--force_forward", action="store_true", default=True,
                        help="RC the sequence when target gene is on minus strand (default: on).")
        sp.add_argument("--no_force_forward", action="store_false", dest="force_forward")

    def add_mitoz(sp: argparse.ArgumentParser) -> None:
        g = sp.add_argument_group("MitoZ")
        g.add_argument("--conda_exe",       default="conda")
        g.add_argument("--conda_env_mitoz", default=None,
                       help="Conda env containing MitoZ (e.g. mitoz3.6).")
        g.add_argument("--mitoz_path",      default=None,
                       help="Path to MitoZ executable or directory.")
        # BUG-V0.29-E: old parser accepted zero/negative values and passed
        # invalid thread counts to MitoZ or ThreadPoolExecutor.
        # g.add_argument("--threads",      type=int, default=16,
        g.add_argument("--threads",      type=positive_int_arg, default=16,
                       help="Threads per MitoZ job (default: 16).")
        # g.add_argument("--max_tasks",    type=int, default=16,
        g.add_argument("--max_tasks",    type=positive_int_arg, default=16,
                       help="Max concurrent annotation jobs (default: 16).")
        g.add_argument("--clade",        default=_CLADE_DEFAULT,
                       help=f"MitoZ --clade value (default: '{_CLADE_DEFAULT}').\n"
                            "⚠ Change this for non-Annelida datasets.")
        g.add_argument("--genetic_code", type=int, default=5,
                       help="Genetic code (default: 5, invertebrate mitochondrial).")
        g.add_argument("--species_name", default=None)
        g.add_argument("--mitoz_extra",  action="append", default=[],
                       metavar="ARG",
                       help="Extra arguments passed verbatim to MitoZ annotate.")

        fq = sp.add_argument_group("FASTQ (optional — for reference-guided annotation)")
        fq.add_argument("--infq", default=None,
                        help="Directory containing paired FASTQ files.")
        fq.add_argument("--fq_position", choices=["flat", "sample_dir"], default="flat")
        fq.add_argument("--suffix_fq", default="_1.clean.fq.gz,_2.clean.fq.gz")

        sch = sp.add_argument_group("Scheduler")
        sch.add_argument("--scheduler_backend",
                         choices=["auto", "snakemake", "threads"], default="auto",
                         help=(
                             "Dispatch backend (default: auto).\n"
                             "  auto      — use Snakemake if available, else threads\n"
                             "  snakemake — always use Snakemake\n"
                             "  threads   — always use ThreadPoolExecutor\n"
                             "MitoZ annotation is parallelised; finalization runs\n"
                             "in the main process regardless of backend."
                         ))
        sch.add_argument("--snakemake_bin",    default=None)
        # BUG-V0.29-E: keep scheduler numeric options positive for predictable
        # Snakemake command generation.
        # sch.add_argument("--smk_cores",        type=int, default=64)
        # sch.add_argument("--smk_mem_mb",       type=int, default=256000)
        # sch.add_argument("--smk_latency_wait", type=int, default=60)
        sch.add_argument("--smk_cores",        type=positive_int_arg, default=64)
        sch.add_argument("--smk_mem_mb",       type=positive_int_arg, default=256000)
        sch.add_argument("--smk_latency_wait", type=positive_int_arg, default=60)
        sch.add_argument("--smk_extra",        action="append", default=[],
                         metavar="FLAG",
                         help="Extra Snakemake flags appended verbatim.\n"
                              "Example: --smk_extra '--profile slurm'")

    # ── annotate ──────────────────────────────────────────────────────────
    a = sub.add_parser("annotate", help="Batch MitoZ annotation only.")
    add_common(a)
    add_mitoz(a)
    a.add_argument("--preserve_phase1_metadata", choices=["yes", "no"], default="no")

    # ── reorient ──────────────────────────────────────────────────────────
    r = sub.add_parser("reorient", help="Reorient GenBank files to start at a target gene.")
    add_common(r)

    # ── run ───────────────────────────────────────────────────────────────
    run_p = sub.add_parser(
        "run",
        help="Full pipeline: annotate → reorient → re-annotate → guaranteed final reorientation.")
    add_common(run_p)
    add_mitoz(run_p)
    run_p.add_argument("--if_reannotation_2", choices=["yes", "no"], default="yes",
                       help="Run Phase-2 re-annotation on reoriented sequences (default: yes).")
    run_p.add_argument("--preserve_phase1_metadata", choices=["yes", "no"], default="yes",
                       help="Restore original GenBank header metadata into final files "
                            "(default: yes; requires GenBank input).")
    run_p.add_argument("--smart_metadata",       choices=["yes", "no"], default="no",
                       help="Enable fuzzy/interactive accession matching (default: no).")
    run_p.add_argument("--auto_accept_metadata", choices=["yes", "no"], default="no",
                       help="Auto-accept best fuzzy match (requires --smart_metadata yes).")

    # ── transfer_metadata ─────────────────────────────────────────────────
    tm = sub.add_parser("transfer_metadata",
                        help="Standalone GenBank header transfer between two files/dirs.")
    tm.add_argument("-m", "--metadata", required=True,
                    help="Original GenBank file/dir (metadata source).")
    tm.add_argument("-n", "--new_body", required=True,
                    help="Re-annotated GenBank file/dir (body source).")
    tm.add_argument("-o", "--out_root", default="transfer_metadata_out")
    tm.add_argument("--if_replace_metadata",  choices=["yes", "no"], default="yes")
    tm.add_argument("--overwrite_new_body",   choices=["yes", "no"], default="no")
    tm.add_argument("--output_name", default=None,
                    help="Output filename (single-file mode only).")
    tm.add_argument("--smart_metadata",       choices=["yes", "no"], default="no")
    tm.add_argument("--auto_accept_metadata", choices=["yes", "no"], default="no")

    # ── replace_organism  (v0.25 / v0.26 / v0.27) ────────────────────────
    ro = sub.add_parser(
        "replace_organism",
        help=(
            "Restore /organism= and SOURCE text in FEATURES source blocks.\n"
            "Fixes MitoZ placeholder values (e.g. 'Test sp.') that survive\n"
            "transfer_metadata because they live in the body, not the header.\n"
            "In 'run' mode this is applied automatically as Step 3.5."
        ),
    )
    ro.add_argument("-m", "--metadata", required=True,
                    help="Original GenBank file/dir (organism/source name source).")
    ro.add_argument("-n", "--new_body", required=True,
                    help="Re-annotated GenBank file/dir whose /organism= needs fixing.")
    ro.add_argument("-o", "--out_root", default="replace_organism_out",
                    help="Output root directory (default: replace_organism_out).")
    ro.add_argument("--overwrite_new_body", choices=["yes", "no"], default="no",
                    help="Overwrite the input files in-place (default: no).")

    return p


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────
def main() -> int:
    if SeqIO is None:
        print("ERROR: Biopython not installed.  pip install biopython", file=sys.stderr)
        return 2

    args = build_parser().parse_args()

    if getattr(args, "conda_env_mitoz", None) == "":
        args.conda_env_mitoz = None

    if getattr(args, "suffix_fq", None):
        parse_suffix_fq(args.suffix_fq)

    if (getattr(args, "auto_accept_metadata", "no") == "yes"
            and getattr(args, "smart_metadata", "no") != "yes"):
        print("ERROR: --auto_accept_metadata yes requires --smart_metadata yes", file=sys.stderr)
        return 2

    if (getattr(args, "clade", None) == _CLADE_DEFAULT
            and args.cmd in ("annotate", "run")):
        print(
            f"\n{'='*68}\n"
            f"  [CLADE_DEFAULT] WARNING\n"
            f"  --clade not set; using default: '{_CLADE_DEFAULT}'\n"
            f"  Set --clade explicitly for non-Annelida datasets.\n"
            f"  See: python3 batch_mitoz.py annotate --help\n"
            f"{'='*68}\n",
            file=sys.stderr)

    dispatch = {
        "annotate":          workflow_annotate,
        "reorient":          workflow_reorient,
        "run":               workflow_run,
        "transfer_metadata": workflow_transfer_metadata,
        "replace_organism":  workflow_replace_organism,
    }
    fn = dispatch.get(args.cmd)
    if fn is None:
        print(f"Unknown command: {args.cmd}", file=sys.stderr)
        return 2
    return fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
