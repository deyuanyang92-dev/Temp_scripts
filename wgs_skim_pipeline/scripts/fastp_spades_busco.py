#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fastp_spades_busco_v6_25.py
===========================
Batch pipeline: fastp (optional) -> SPAdes (auto-resume) -> BUSCO
Extra subcommands: clean | check | status

NEW in v6.17 — Snakemake SPAdes continue-order fix + shell hardening
===============================================
BUG-FIX 1 (CRITICAL — rc=127): Removed '--' from `conda run` prefix in
  generated Snakefile.  Older conda versions (< 23.x) use
  `nargs=REMAINDER` in `conda run` argparse, which does NOT consume '--'
  as an option-terminator; instead it passes '--' as the first token of
  the command to execute.  Bash then searches PATH for a binary literally
  named '--' → "command not found" → rc=127 on every SPAdes/BUSCO job.

BUG-FIX 2 (CRITICAL — empty logs): Removed `mv -f` of stdout_spades.log
  and stderr_spades.log in the spades_wipe_inline shell block.  Snakemake
  opens {log.out} and {log.err} file descriptors BEFORE the shell block
  runs.  Renaming the files under Snakemake's open fds causes Snakemake
  to read empty/stale logs and misreport job status.  The subsequent
  `> {log.out}` redirect truncates the file in-place anyway, so no data
  is lost.  Same fix applied to busco_wipe_inline.

BUG-FIX 3 (output dir scanned as sample): When --input_dir and --out_root
  share a parent (e.g. -i /WGS -o /WGS/results), the output directory
  "results" was picked up as a sample with missing reads.  Now
  _collect_tasks_from_dir() accepts an exclude_dirs parameter and
  collect_samples_from_dir() passes the resolved --out_root to it.

BUG-FIX 4 (relative BUSCO lineage path): Relative paths like
  '../../data/metazoa_odb12/' were passed to the Snakefile verbatim.
  Snakemake runs from --out_root as CWD, so relative paths resolve
  incorrectly → BUSCO "lineage not found".  Now resolved to absolute
  before Snakefile generation.  Same fix for --busco_fasta_dir.

BUG-FIX 5 (double status table): run_snakemake_pipeline() printed the
  status table a second time after main() already did.  Removed duplicate.

NEW in v6.24 — Gemini/Claude joint audit fixes
===============================================
BUG-FIX 6 (CRITICAL — Snakemake SyntaxError): rule busco was emitting an
  empty "input:" block when no SPAdes dependency existed (busco-standalone
  mode).  Snakemake treats an empty input: as a syntax error.  Fixed by
  emitting the "input:" keyword only when busco_spades_dep is non-empty.

BUG-FIX 7 (CRITICAL — dead code): SPAdes --continue was never invoked in
  the Snakemake backend.  The v6.10 template refactor stored the fresh
  command in a Bash variable (FRESH_CMD) and always ran it unconditionally,
  silently discarding the checkpoint-detect → --continue → fallback-fresh
  logic.  Restored to match the ThreadPool backend (run_spades_one).

BUG-FIX 8 (CRITICAL — word splitting): FRESH_CMD="..." stored the SPAdes
  command in a plain Bash string variable.  Arguments containing spaces
  (e.g. --spades_extra "--custom 'my value'") were split by Bash variable
  expansion, passing broken flags to SPAdes.  Fixed by injecting fresh_cmd
  and cont_cmd directly as literal shell text (same pattern as rule fastp
  and rule busco already used).

BUG-FIX 9 (biological): WGS BUSCO default fasta was contigs.fasta.
  Standard practice is scaffolds.fasta (more contiguous, higher N50),
  which gives accurate completeness scores.  Changed AUTO default in
  choose_busco_fasta(), _smk_busco_fasta_name(), and
  _busco_fasta_name_for_display().  Fallback chain already preferred
  scaffolds via spades_expected_outputs() — AUTO default now consistent.

BUG-FIX 10 (portability): rule busco shell used compgen -G to test for
  partial-state run_*/ directories.  compgen is a bash builtin; Snakemake
  defaults to /bin/sh which is dash on Debian/Ubuntu → rc=127 crash.
  Replaced with a POSIX-compatible for-glob + [ -e ] guard loop.

BUG-FIX 11 (silent data error): choose_busco_fasta() silently fell through
  to the SPAdes output directory when --busco_fasta_dir was set but the
  expected fasta was missing.  This could run BUSCO on an unrelated
  assembly without any warning.  Now returns the expected (missing) path
  and emits a WARNING, so the downstream existence check surfaces a clear
  error message.

NEW in v6.25 — check/status subcommand bug-fixes
===============================================
BUG-FIX 12 (check blocked by empty-samples): preflight_env_check() was
  delegating to preflight_full_check() with fake_samples=[]. _chk_samples()
  treats an empty list as a fatal error ("No samples found"), making `check`
  always fail before it could inspect any tools.  Fixed by adding a dedicated
  preflight_tools_only_check() that runs only checks A–E (Python, conda,
  tools, snakemake, resources) with no sample or busco_lineage requirements.

BUG-FIX 13 (check blocked by missing busco_lineage): _chk_params() requires
  --busco_lineage whenever cmd in ("all","busco","spades2busco").  The old
  `check` handler derived check_cmd="all" when all three stages were selected
  (the default), hitting the lineage guard even though the user never intended
  to supply one.  Fixed: `check` now calls preflight_tools_only_check() which
  skips parameter-sanity checks entirely.

BUG-FIX 14 (status auto-scan wrong directory): The status auto-scan
  (no -i / --sample_table) only looked under out_root/spades_dir.  Running
  `spades` or `busco` subcommands stores outputs directly under out_root/<sid>/
  (direct-output mode), so those samples were invisible to `status`.  The
  scan now checks spades_dir, then busco_dir, then falls back to scanning
  out_root directly — covering all output layouts.

BUG-FIX 15 (status sample IDs include suffix): The auto-scan appended
  d.name as the sample ID without stripping --suffix_spades_out (e.g.
  _k77).  Status table IDs then differed from the IDs stored in stamp
  files → every sample falsely showed TODO.  Now strips the suffix.

BUG-FIX 16 (status/busco scan picks up non-sample dirs): Directory scans
  iterated every subdirectory of out_root without filtering, treating
  pipeline-internal dirs (01.fastp, .snakemake, __pycache__, …) as
  sample IDs.  All scans now skip hidden dirs, known pipeline-dir names,
  and other non-sample directories.

NEW in v6.25.1 — status_cmd inference + shim safety
===============================================
BUG-FIX 17 (status overall always wrong): `status` hardcoded
  status_cmd="all", so scan_all_status() always treated fastp/spades/busco
  as all active.  When checking results of a `spades`-only run (no fastp),
  the fastp column showed TODO and dragged "overall" to TODO even for fully
  assembled samples.  Fixed: status_cmd is now auto-inferred from the
  presence of 01.fastp/, 02.spades/, 03.busco/ directories (and direct-
  output stamp files), then mapped to the correct subcommand name.
  A new --pipeline_cmd flag lets users override the inference explicitly.

BUG-FIX 18 (preflight_env_check shim crashes): The deprecated
  preflight_env_check() backward-compat shim still delegated to
  preflight_full_check() with fake_samples=[].  Any external caller
  (e.g. wrapper scripts) would hit the _chk_samples() "No samples found"
  fatal error instantly.  Fixed: shim now delegates to
  preflight_tools_only_check() — same as the `check` subcommand does —
  so it never requires samples or busco_lineage.

All v6.10/v6.7/v6.4/v6.3 features retained."""

from __future__ import annotations

import argparse
import json as _json
import logging
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import textwrap
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


# ==============================================================
# 0. Python-version compatibility shim
# ==============================================================

if not hasattr(argparse, "BooleanOptionalAction"):
    class _BooleanOptionalAction(argparse.Action):
        """Backport of argparse.BooleanOptionalAction for Python < 3.9."""

        def __init__(self, option_strings, dest, default=None,
                     type=None, choices=None, required=False,
                     help=None, metavar=None):
            _opts: List[str] = []
            self._pos_opts: Set[str] = set()
            self._neg_opts: Set[str] = set()
            for opt in option_strings:
                _opts.append(opt)
                self._pos_opts.add(opt)
                if opt.startswith("--"):
                    neg = "--no-" + opt[2:]
                    _opts.append(neg)
                    self._neg_opts.add(neg)
                    neg2 = "--no_" + opt[2:]
                    if neg2 != neg:
                        _opts.append(neg2)
                        self._neg_opts.add(neg2)
            super().__init__(
                option_strings=_opts, dest=dest, nargs=0,
                default=default, type=type, choices=choices,
                required=required, help=help, metavar=metavar,
            )

        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest,
                    option_string not in self._neg_opts)

        def format_usage(self):
            return " | ".join(self.option_strings)

    argparse.BooleanOptionalAction = _BooleanOptionalAction


# ==============================================================
# 1. Logging
# ==============================================================

def setup_root_logger(log_file: Optional[Path] = None,
                      level: int = logging.INFO) -> None:
    root = logging.getLogger()
    root.setLevel(level)
    for h in list(root.handlers):
        root.removeHandler(h)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    root.addHandler(sh)
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(str(log_file), encoding="utf-8")
        fh.setFormatter(fmt)
        root.addHandler(fh)


# ==============================================================
# 2. General helpers
# ==============================================================

def safe_name(s: str) -> str:
    s = str(s).strip().replace(" ", "_")
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", s)
    return s.strip("_") or "NA"


def which_or_none(exe: str) -> Optional[str]:
    try:
        return shutil.which(exe)
    except Exception:
        return None


def tail_text(s: str, n: int = 120) -> str:
    lines = (s or "").splitlines()
    return "\n".join(lines[-n:])


def extra_args_from_list(lst: Optional[List[str]]) -> List[str]:
    out: List[str] = []
    for s in lst or []:
        if s is not None:
            out.extend(shlex.split(str(s)))
    return out


def resolve_list_arg(raw: Optional[List[str]],
                     base_dir: Path) -> Optional[List[str]]:
    if not raw:
        return None
    if len(raw) == 1:
        p = Path(raw[0])
        if not p.is_absolute():
            p2 = base_dir / raw[0]
            if p2.is_file():
                p = p2
        if p.is_file():
            items: List[str] = []
            seen: Set[str] = set()
            for line in p.read_text(encoding="utf-8",
                                     errors="ignore").splitlines():
                s = line.strip()
                if s and not s.startswith("#") and s not in seen:
                    seen.add(s)
                    items.append(s)
            return items
    items = []
    seen_raw: Set[str] = set()
    for v in raw:
        for x in str(v).split(","):
            x = x.strip()
            if x and x not in seen_raw:
                seen_raw.add(x)
                items.append(x)
    return items or None


def resolve_exe_from_path(exe_name: str,
                          tool_path: Optional[str]) -> Optional[Path]:
    if not tool_path:
        return None
    p = Path(tool_path).expanduser()
    return (p / exe_name).resolve() if p.is_dir() else p.resolve()


def _colocated_python(script_path: Path) -> str:
    script_bin = script_path.parent
    for py_name in ("python3", "python"):
        co_py = script_bin / py_name
        if co_py.exists() and co_py.is_file():
            return str(co_py)
    return sys.executable


def build_tool_cmd(*, tool_label: str, exe_name: str,
                   base_args: List[str],
                   conda_env: Optional[str],
                   tool_path: Optional[str],
                   conda_exe: str = "conda") -> List[str]:
    if conda_env:
        conda_bin = which_or_none(conda_exe)
        if not conda_bin:
            raise FileNotFoundError(
                f"[{tool_label}] '{conda_exe}' not found in PATH.\n"
                f"  Hint: install conda/mamba, or specify --conda_exe /path/to/conda")
        # IMPORTANT: do NOT add '--' after `conda run`.
        # Older conda versions can pass that literal token through to bash,
        # producing: `--: command not found`.
        return ([conda_bin, "run", "-n", conda_env,
                 "--no-capture-output", exe_name] + base_args)
    exe_path = resolve_exe_from_path(exe_name, tool_path)
    if exe_path is not None:
        if not exe_path.exists():
            raise FileNotFoundError(
                f"[{tool_label}] binary not found at: {exe_path}\n"
                f"  Hint: check --{tool_label.lower()}_path value")
        if exe_path.suffix.lower() == ".py":
            return [_colocated_python(exe_path), str(exe_path)] + base_args
        return [str(exe_path)] + base_args
    found = which_or_none(exe_name)
    if not found:
        raise FileNotFoundError(
            f"[{tool_label}] '{exe_name}' not found in PATH.\n"
            f"  Fix options (pick one):\n"
            f"    1. conda activate <env_with_{exe_name}>\n"
            f"    2. --conda_env_{tool_label.lower()} ENV_NAME\n"
            f"    3. --{tool_label.lower()}_path /path/to/{exe_name}\n"
            f"    4. export PATH=/dir/containing/{exe_name}:$PATH")
    found_path = Path(found)
    if found_path.suffix.lower() == ".py":
        return [_colocated_python(found_path), found] + base_args
    return [found] + base_args


# -- stamp helpers -----------------------------------------------

def write_done(stamp: Path) -> None:
    stamp.parent.mkdir(parents=True, exist_ok=True)
    stamp.write_text("ok=1\nreason=ok\n", encoding="utf-8")


def done_ok(stamp: Path) -> bool:
    if not stamp.exists():
        return False
    try:
        return any(ln.strip() == "ok=1"
                   for ln in stamp.read_text(
                       encoding="utf-8", errors="ignore").splitlines()[:3])
    except Exception:
        return False


def write_fail(fail_path: Path, reason: str) -> None:
    fail_path.parent.mkdir(parents=True, exist_ok=True)
    fail_path.write_text(f"ok=0\nreason={reason}\n", encoding="utf-8")


def _run_cmd(cmd: List[str], cwd: Optional[Path] = None,
             dryrun: bool = False) -> Tuple[int, str, str]:
    if dryrun:
        logging.info(f"[DRYRUN] {' '.join(str(c) for c in cmd)}")
        return 0, "", ""
    # --- v6: transparent command logging ---
    cmd_str = " ".join(str(c) for c in cmd)
    logging.info(f"[CMD] {cmd_str}")
    try:
        p = subprocess.run(cmd, cwd=str(cwd) if cwd else None,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           text=True)
        return int(p.returncode), p.stdout or "", p.stderr or ""
    except FileNotFoundError as e:
        return 127, "", f"Command not found: {e}"
    except PermissionError as e:
        return 126, "", f"Permission denied: {e}"
    except Exception as e:
        return 1, "", f"Unexpected error: {e}"


# ==============================================================
# 2a. System resource detection + cap
# ==============================================================

def get_system_resources() -> Tuple[int, int]:
    cpus = os.cpu_count() or 1
    mem_gb = 32
    # Linux: /proc/meminfo
    try:
        with open("/proc/meminfo", encoding="utf-8") as f:
            for line in f:
                if line.startswith("MemTotal"):
                    mem_kb = int(line.split()[1])
                    mem_gb = max(1, mem_kb // (1024 * 1024))
                    return cpus, mem_gb   # found — return immediately
    except Exception:
        pass
    # macOS: sysctl hw.memsize
    try:
        r = subprocess.run(
            ["sysctl", "-n", "hw.memsize"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, timeout=5)
        if r.returncode == 0:
            mem_bytes = int(r.stdout.strip())
            mem_gb = max(1, mem_bytes // (1024 ** 3))
    except Exception:
        pass
    return cpus, mem_gb


def get_disk_free_gb(path: Path) -> float:
    try:
        st = os.statvfs(str(path))
        return (st.f_bavail * st.f_frsize) / (1024 ** 3)
    except Exception:
        return -1.0


def cap_resources(params: dict, max_pct: float = 0.90) -> dict:
    cpus, mem_gb = get_system_resources()
    max_threads = max(1, int(cpus * max_pct))
    max_mem = max(1, int(mem_gb * max_pct))
    params = dict(params)
    for key in ("threads_fastp", "threads_spades", "threads_busco"):
        requested = int(params.get(key, 0) or 0)
        if requested > max_threads:
            logging.warning(
                f"[RESOURCE_CAP] {key}: {requested} -> {max_threads} "
                f"(<= {int(max_pct * 100)}% of {cpus} CPUs)")
            params[key] = max_threads
    mem_key = "spades_mem_gb"
    requested_mem = int(params.get(mem_key, 0) or 0)
    if requested_mem > max_mem:
        logging.warning(
            f"[RESOURCE_CAP] {mem_key}: {requested_mem} -> {max_mem} GB "
            f"(<= {int(max_pct * 100)}% of {mem_gb} GB RAM)")
        params[mem_key] = max_mem
    par = int(params.get("max_parallel", 8) or 8)
    if par > max_threads:
        logging.warning(f"[RESOURCE_CAP] max_parallel: {par} -> {max_threads}")
        params["max_parallel"] = max_threads
    logging.info(
        f"[RESOURCE] system  CPUs={cpus}  RAM={mem_gb} GB  "
        f"cap={int(max_pct * 100)}%  -> "
        f"max_threads={max_threads}  max_mem_gb={max_mem}")
    return params


# ==============================================================
# 2b. Per-sample directory helpers
# ==============================================================

def spades_sample_dir(out_root: Path, params: dict, sid: str) -> Path:
    suffix = params.get("suffix_spades_out", "") or ""
    if params.get("direct_spades_output", False):
        return out_root / f"{sid}{suffix}"
    return out_root / params["spades_dir"] / f"{sid}{suffix}"


def busco_run_dir(out_root: Path, params: dict, sid: str) -> Path:
    if params.get("direct_busco_output", False):
        return out_root / sid
    return out_root / params["busco_dir"] / sid


# ==============================================================
# 2c. fastp completeness verification
# ==============================================================

def fastp_output_complete(stage_dir: Path, sid: str,
                          done_stamp: Path) -> bool:
    if not done_ok(done_stamp):
        return False
    r1 = stage_dir / f"{sid}.R1.clean.fastq.gz"
    r2 = stage_dir / f"{sid}.R2.clean.fastq.gz"
    if not r1.exists() or r1.stat().st_size == 0:
        return False
    if not r2.exists() or r2.stat().st_size == 0:
        return False
    json_p = stage_dir / "fastp.json"
    if not json_p.exists():
        return False
    try:
        data = _json.loads(json_p.read_text(encoding="utf-8", errors="ignore"))
        return "summary" in data
    except Exception:
        return False


def _invalidate_fastp_output(stage_dir: Path, sid: str,
                              done_stamp: Path, fail_stamp: Path,
                              reason: str) -> None:
    done_stamp.unlink(missing_ok=True)
    for fname in [f"{sid}.R1.clean.fastq.gz",
                  f"{sid}.R2.clean.fastq.gz",
                  "fastp.json", "fastp.html"]:
        p = stage_dir / fname
        if p.exists():
            try:
                p.unlink()
            except Exception:
                pass
    write_fail(fail_stamp, reason)
    logging.warning(f"[FASTP_CORRUPT] {sid}: {reason} -- removed partial output")


# ==============================================================
# 2d. BUSCO version detection + partial-state helpers
# ==============================================================

_BUSCO_VERSION_CACHE: Dict[str, Tuple[int, int]] = {}


def detect_busco_version(busco_exe: str,
                          conda_env: Optional[str],
                          tool_path: Optional[str],
                          conda_exe: str) -> Tuple[int, int]:
    cache_key = f"{conda_env or ''}::{tool_path or ''}::{busco_exe}"
    if cache_key in _BUSCO_VERSION_CACHE:
        return _BUSCO_VERSION_CACHE[cache_key]
    ver: Tuple[int, int] = (5, 0)
    try:
        cmd = build_tool_cmd(
            tool_label="BUSCO", exe_name=busco_exe,
            base_args=["--version"],
            conda_env=conda_env, tool_path=tool_path, conda_exe=conda_exe)
        p = subprocess.run(cmd, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, text=True, timeout=30)
        txt = (p.stdout + p.stderr).strip()
        m = re.search(r"(\d+)\.(\d+)", txt)
        if m:
            ver = (int(m.group(1)), int(m.group(2)))
            logging.info(f"[BUSCO_VER] detected BUSCO {ver[0]}.{ver[1]}")
    except Exception as exc:
        logging.debug(f"[BUSCO_VER] version detection failed: {exc}")
    _BUSCO_VERSION_CACHE[cache_key] = ver
    return ver


def busco_has_partial_state(run_dir: Path) -> bool:
    """Return True when BUSCO has a real checkpoint to restart from.

    Official BUSCO docs on -r/--restart:
      "Continue a run that had already partially completed.
       Restarting skips calls to tools that have completed but
       performs all pre- and post-processing steps."

    A real restartable checkpoint requires run_<lineage>/ to exist
    and contain actual analysis results (not just empty dirs or logs).
    The run_*/ directory is written by BUSCO during analysis.
    An empty run_*/ or one with only log files is NOT a real checkpoint.
    """
    if not run_dir.is_dir():
        return False
    if busco_has_summary(run_dir):
        return False   # already complete, not partial
    # Check for non-empty run_*/ directory = real checkpoint
    for run_d in run_dir.glob("run_*/"):
        if run_d.is_dir():
            try:
                children = list(run_d.iterdir())
                # Must have more than just log files
                non_log = [c for c in children
                           if c.name not in ("logs", "busco_db")
                           and not c.name.endswith(".log")]
                if non_log:
                    return True
            except PermissionError:
                pass
    return False


def _busco_build_restart_args(ver: Tuple[int, int],
                               run_dir: Path) -> List[str]:
    return ["--restart"]


# ==============================================================
# 2e. Skip-set loader
# ==============================================================

def load_skip_set(raw_vals: Optional[List[str]],
                  base_dir: Path,
                  summary_tsv: Optional[Path] = None,
                  stage: str = "all") -> Set[str]:
    skip: Set[str] = set()
    ids = resolve_list_arg(raw_vals, base_dir) if raw_vals else None
    if ids:
        skip.update(safe_name(x) for x in ids)
    if summary_tsv and summary_tsv.is_file():
        col_map: Dict[str, int] = {}
        try:
            for line in summary_tsv.read_text(
                    encoding="utf-8", errors="ignore").splitlines():
                stripped = line.strip()
                if not stripped:
                    continue
                parts = stripped.split("\t")
                if not col_map:
                    for ci, h in enumerate(parts):
                        col_map[h.lower().lstrip("#").strip()] = ci
                    continue
                sid_ci = col_map.get("sample", 0)
                sid_val = safe_name(parts[sid_ci]) if sid_ci < len(parts) else ""
                if not sid_val:
                    continue

                def _done_col(name: str) -> bool:
                    ci = col_map.get(name)
                    return (ci is not None and ci < len(parts)
                            and parts[ci].strip() == "1")

                if stage == "fastp" and _done_col("fastp_done"):
                    skip.add(sid_val)
                elif stage == "spades" and _done_col("spades_done"):
                    skip.add(sid_val)
                elif stage == "busco" and _done_col("busco_done"):
                    skip.add(sid_val)
                elif stage == "all" and all(
                    _done_col(c) for c in
                    ["fastp_done", "spades_done", "busco_done"]
                ):
                    skip.add(sid_val)
        except Exception as exc:
            logging.warning(f"[SKIP_SET] could not parse {summary_tsv}: {exc}")
    return skip


# ==============================================================
# 2f. Pre-run result scanners  (v6.1: backward-compatible with v5)
# ==============================================================

def _spades_is_complete(sd: Path, data_type: str) -> bool:
    """Check if SPAdes output is complete -- works with or without .done stamp."""
    if done_ok(sd / ".spades.done"):
        return True
    # v5 backward compat: check raw output files even without stamp
    return bool(spades_success_files(sd, data_type))


def _busco_is_complete(rd: Path) -> bool:
    """Check if BUSCO output is complete -- works with or without .done stamp."""
    if busco_has_summary(rd) and done_ok(rd / ".busco.done"):
        return True
    # v5 backward compat: check summary file even without stamp
    if busco_has_summary(rd):
        return True
    return False


def _fastp_is_complete(stage_dir: Path, sid: str) -> bool:
    """Check if fastp output is complete -- works with or without stamp."""
    done_stamp = stage_dir / ".fastp.done"
    if fastp_output_complete(stage_dir, sid, done_stamp):
        return True
    # v5 backward compat: check output files even without stamp
    r1 = stage_dir / f"{sid}.R1.clean.fastq.gz"
    r2 = stage_dir / f"{sid}.R2.clean.fastq.gz"
    json_p = stage_dir / "fastp.json"
    if r1.exists() and r1.stat().st_size > 0 and \
       r2.exists() and r2.stat().st_size > 0 and json_p.exists():
        try:
            data = _json.loads(json_p.read_text(encoding="utf-8", errors="ignore"))
            if "summary" in data:
                # Retroactively write stamp for future runs
                write_done(done_stamp)
                return True
        except Exception:
            pass
    return False


def prescan_spades(samples: List["Sample"], out_root: Path,
                   params: dict) -> Tuple[List["Sample"], List["Sample"]]:
    done_list: List["Sample"] = []
    todo_list: List["Sample"] = []
    data_type = params.get("data_type", "wgs")
    for s in samples:
        sd = spades_sample_dir(out_root, params, s.sid)
        if _spades_is_complete(sd, data_type):
            done_list.append(s)
        else:
            todo_list.append(s)
    if done_list:
        sids = ", ".join(s.sid for s in done_list[:10])
        logging.info(
            f"[PRESCAN_SPADES] {len(done_list)} already assembled -- skip: "
            f"{sids}{'...' if len(done_list) > 10 else ''}")
    logging.info(
        f"[PRESCAN_SPADES] done={len(done_list)}  to_run={len(todo_list)}  "
        f"total={len(samples)}")
    return done_list, todo_list


def prescan_busco(samples: List["Sample"], out_root: Path,
                  params: dict) -> Tuple[List["Sample"], List["Sample"]]:
    done_list: List["Sample"] = []
    todo_list: List["Sample"] = []
    for s in samples:
        rd = busco_run_dir(out_root, params, s.sid)
        if _busco_is_complete(rd):
            done_list.append(s)
        else:
            todo_list.append(s)
    if done_list:
        sids = ", ".join(s.sid for s in done_list[:10])
        logging.info(
            f"[PRESCAN_BUSCO] {len(done_list)} already done -- skip: "
            f"{sids}{'...' if len(done_list) > 10 else ''}")
    logging.info(
        f"[PRESCAN_BUSCO] done={len(done_list)}  to_run={len(todo_list)}  "
        f"total={len(samples)}")
    return done_list, todo_list


def prescan_fastp(samples: List["Sample"], out_root: Path,
                  params: dict) -> Tuple[List["Sample"], List["Sample"]]:
    done_list: List["Sample"] = []
    todo_list: List["Sample"] = []
    for s in samples:
        stage_dir = out_root / params["fastp_dir"] / s.sid
        if _fastp_is_complete(stage_dir, s.sid):
            done_list.append(s)
        else:
            todo_list.append(s)
    if done_list:
        sids = ", ".join(s.sid for s in done_list[:10])
        logging.info(
            f"[PRESCAN_FASTP] {len(done_list)} already cleaned -- skip: "
            f"{sids}{'...' if len(done_list) > 10 else ''}")
    logging.info(
        f"[PRESCAN_FASTP] done={len(done_list)}  to_run={len(todo_list)}  "
        f"total={len(samples)}")
    return done_list, todo_list


# ==============================================================
# 2f-2. Full status scan & status table  (NEW v6.5)
# ==============================================================

# Status categories for each sample×stage
_ST_DONE    = "DONE"      # complete
_ST_PARTIAL = "PARTIAL"   # partial / interrupted — highest priority
_ST_TODO    = "TODO"      # not started
_ST_SKIP    = "SKIP"      # stage not requested for this run


def _stage_status(stage: str, params: dict,
                  out_root: Path, s: "Sample") -> str:
    """Return one of _ST_* for *stage* of sample *s*."""
    data_type = params.get("data_type", "wgs")
    if stage == "fastp":
        sd = out_root / params["fastp_dir"] / s.sid
        if _fastp_is_complete(sd, s.sid):
            return _ST_DONE
        r1 = sd / f"{s.sid}.R1.clean.fastq.gz"
        r2 = sd / f"{s.sid}.R2.clean.fastq.gz"
        if r1.exists() or r2.exists():
            return _ST_PARTIAL
        return _ST_TODO

    if stage == "spades":
        sd = spades_sample_dir(out_root, params, s.sid)
        if _spades_is_complete(sd, data_type):
            return _ST_DONE
        if spades_has_partial_state(sd):
            return _ST_PARTIAL
        return _ST_TODO

    if stage == "busco":
        rd = busco_run_dir(out_root, params, s.sid)
        if _busco_is_complete(rd):
            return _ST_DONE
        if busco_has_partial_state(rd):
            return _ST_PARTIAL
        return _ST_TODO

    return _ST_SKIP


def scan_all_status(
    samples: List["Sample"],
    out_root: Path,
    params: dict,
    cmd: str,
) -> List[Dict[str, str]]:
    """
    Scan every sample for every stage and return a list of dicts:
      {sid, fastp, spades, busco, overall, priority}
    overall: DONE | RESUME | TODO
    priority: "10" (resume) | "0" (todo) | "-" (done, excluded from run)
    """
    run_fastp  = cmd in ("all", "fastp")  and params.get("use_fastp", True)
    run_spades = cmd in ("all", "spades", "spades2busco")
    run_busco  = cmd in ("all", "busco",  "spades2busco")

    rows: List[Dict[str, str]] = []
    for s in samples:
        fs = _stage_status("fastp",  params, out_root, s) if run_fastp  else _ST_SKIP
        ss = _stage_status("spades", params, out_root, s) if run_spades else _ST_SKIP
        bs = _stage_status("busco",  params, out_root, s) if run_busco  else _ST_SKIP

        active = [v for v in (fs, ss, bs) if v != _ST_SKIP]
        if all(v == _ST_DONE for v in active):
            overall = _ST_DONE
            prio    = "-"
        elif any(v == _ST_PARTIAL for v in active):
            overall = "RESUME"
            prio    = "10"
        else:
            overall = _ST_TODO
            prio    = "0"

        rows.append({"sid": s.sid,
                     "fastp":  fs,
                     "spades": ss,
                     "busco":  bs,
                     "overall": overall,
                     "priority": prio})
    return rows


def print_status_table(rows: List[Dict[str, str]],
                       cmd: str, params: dict,
                       out_root: Path,
                       max_rows: int = 60) -> None:
    """Print a compact status table to stdout."""
    run_fastp  = cmd in ("all", "fastp")  and params.get("use_fastp", True)
    run_spades = cmd in ("all", "spades", "spades2busco")
    run_busco  = cmd in ("all", "busco",  "spades2busco")

    # Colour / symbol map (ANSI, degrades gracefully in plain terminals)
    _SYM = {
        _ST_DONE:    "✔ DONE   ",
        _ST_PARTIAL: "↻ PARTIAL",
        _ST_TODO:    "· TODO   ",
        _ST_SKIP:    "  —      ",
    }

    total   = len(rows)
    n_done  = sum(1 for r in rows if r["overall"] == _ST_DONE)
    n_res   = sum(1 for r in rows if r["overall"] == "RESUME")
    n_todo  = sum(1 for r in rows if r["overall"] == _ST_TODO)

    sep = "─" * 72
    print("\n" + "═"*72)
    print(f"  TASK STATUS SCAN   output: {out_root}   total samples: {total}")
    print("═"*72)
    print(f"  {'Sample':<32s}  {'fastp':<9}  {'spades':<9}  {'busco':<9}  overall")
    print(sep)

    # Sort: RESUME first, then TODO, then DONE
    _order = {"RESUME": 0, _ST_TODO: 1, _ST_DONE: 2}
    sorted_rows = sorted(rows, key=lambda r: _order.get(r["overall"], 9))

    printed = 0
    for r in sorted_rows:
        sid_disp = r["sid"] if len(r["sid"]) <= 32 else r["sid"][:29] + "..."
        f_sym = _SYM[r["fastp"]]  if run_fastp  else "  —      "
        s_sym = _SYM[r["spades"]] if run_spades else "  —      "
        b_sym = _SYM[r["busco"]]  if run_busco  else "  —      "
        print(f"  {sid_disp:<32s}  {f_sym}  {s_sym}  {b_sym}  {r['overall']}")
        printed += 1
        if printed >= max_rows and total > max_rows:
            print(f"  ... (showing {max_rows}/{total}; use `status` subcommand for full list)")
            break

    print(sep)
    print(f"  Summary:  ✔ DONE={n_done}  ↻ RESUME(priority)={n_res}  · TODO={n_todo}")
    print(f"  Will run: {n_res + n_todo} samples  (skip {n_done} already complete)")
    print("═"*72 + "\n")


# ==============================================================
# 2g. BUSCO local-lineage detection (NEW v6)
# ==============================================================

def is_local_busco_lineage(lineage: str) -> bool:
    """Return True when *lineage* is a local path (directory or tarball)."""
    if not lineage:
        return False
    p = Path(lineage)
    if p.is_dir():
        return True
    if p.is_file() and p.suffix in (".gz", ".tar", ".tgz"):
        return True
    if os.sep in lineage or "/" in lineage:
        return True
    return False


# ==============================================================
# 2h. Execution plan + user confirmation  (NEW v6.1)
# ==============================================================

def _busco_fasta_name_for_display(params: dict) -> str:
    """Return the fasta filename that BUSCO will consume (for display)."""
    data_type = params.get("data_type", "wgs")
    which = (params.get("which_fasta_to_busco") or "").strip()
    if not which or which.upper() == "AUTO":
        return "transcripts.fasta" if data_type == "rna" else "scaffolds.fasta"
    return which


def print_execution_plan(cmd: str, params: dict,
                         samples: List["Sample"],
                         out_root: Path,
                         use_smk: bool) -> dict:
    """Print a detailed execution plan. Returns prescan counts dict."""
    cpus, mem_gb = get_system_resources()
    free_gb = get_disk_free_gb(out_root)
    data_type = params.get("data_type", "wgs")
    total = len(samples)

    # ---- prescan existing results ----
    fastp_done = spades_done = busco_done = 0
    fastp_todo = spades_todo = busco_todo = total

    if cmd in ("all", "fastp") and params.get("use_fastp", True):
        for s in samples:
            sd = out_root / params["fastp_dir"] / s.sid
            if _fastp_is_complete(sd, s.sid):
                fastp_done += 1
        fastp_todo = total - fastp_done

    if cmd in ("all", "spades", "spades2busco"):
        for s in samples:
            sd = spades_sample_dir(out_root, params, s.sid)
            if _spades_is_complete(sd, data_type):
                spades_done += 1
        spades_todo = total - spades_done

    if cmd in ("all", "busco", "spades2busco"):
        for s in samples:
            rd = busco_run_dir(out_root, params, s.sid)
            if _busco_is_complete(rd):
                busco_done += 1
        busco_todo = total - busco_done

    # ---- BUSCO fasta preview ----
    busco_fasta_display = _busco_fasta_name_for_display(params)
    if params.get("busco_fasta_dir"):
        busco_fasta_source = f"from --busco_fasta_dir ({params['busco_fasta_dir']})"
    else:
        busco_fasta_source = f"from SPAdes output dir"

    # ---- print ----
    sep = "=" * 60
    print(f"\n{sep}")
    print(f"  PIPELINE EXECUTION PLAN  (v6.25)")
    print(f"{sep}")
    print()

    # -- basic info --
    print(f"  BASIC INFORMATION")
    print(f"    Subcommand       : {cmd}")
    print(f"    Total Samples    : {total}")
    print(f"    Data Type        : {data_type}")
    print(f"    Output Root      : {out_root}")
    print(f"    Backend          : {'Snakemake' if use_smk else 'ThreadPool'}")
    print()

    # -- pipeline stages --
    print(f"  PIPELINE STAGES")
    if cmd in ("all", "fastp") and params.get("use_fastp", True):
        print(f"    [1] FASTP  (Quality Control)")
    if cmd in ("all", "spades", "spades2busco"):
        print(f"    [2] SPAdes (Assembly)")
    if cmd in ("all", "busco", "spades2busco"):
        print(f"    [3] BUSCO  (Completeness)")
    print()

    # -- key configurations --
    print(f"  KEY CONFIGURATIONS")
    if cmd in ("all", "fastp") and params.get("use_fastp", True):
        print(f"    FASTP:")
        print(f"      threads          : {params['threads_fastp']}")
        print(f"      adapter_detect   : {params.get('fastp_detect_adapter_for_pe', True)}")
        extra = extra_args_from_list(params.get("fastp_extra", []))
        if extra:
            print(f"      extra args       : {' '.join(extra)}")
        print(f"      [adjust] --threads_fastp N / --fastp_extra '...'")
        print()

    if cmd in ("all", "spades", "spades2busco"):
        asm = params.get("assembly_mode") or "(default)"
        print(f"    SPAdes:")
        print(f"      threads          : {params['threads_spades']}")
        print(f"      memory           : {params['spades_mem_gb']} GB")
        print(f"      assembly_mode    : {asm}")
        print(f"      auto_clean       : {params.get('spades_clean', True)}")
        sfx = params.get("suffix_spades_out", "")
        if sfx:
            print(f"      suffix_out       : {sfx}")
        extra = extra_args_from_list(params.get("spades_extra", []))
        if extra:
            print(f"      extra args       : {' '.join(extra)}")
        print(f"      [adjust] --threads_spades N / --spades_mem_gb N / --assembly_mode X")
        print()

    if cmd in ("all", "busco", "spades2busco"):
        print(f"    BUSCO:")
        print(f"      lineage          : {params.get('busco_lineage', '')}")
        print(f"      mode             : {params.get('busco_mode') or ('transcriptome' if data_type == 'rna' else 'genome')}")
        print(f"      threads          : {params['threads_busco']}")
        print(f"      offline          : {params.get('busco_offline', False)}")
        print(f"      input fasta      : {busco_fasta_display}  ({busco_fasta_source})")
        extra = extra_args_from_list(params.get("busco_extra", []))
        if extra:
            print(f"      extra args       : {' '.join(extra)}")
        print(f"      [adjust] --busco_lineage X / --which_fasta_to_busco X / --busco_mode X")
        print()

    # -- resource allocation --
    print(f"  SYSTEM RESOURCES")
    print(f"    CPUs available   : {cpus}")
    print(f"    RAM available    : {mem_gb} GB")
    if free_gb >= 0:
        print(f"    Disk free        : {free_gb:.1f} GB")
    print(f"    max_parallel     : {params['max_parallel']}")
    print(f"    resource_cap     : {params.get('resource_cap_pct', 90)}%")
    print(f"    [adjust] --max_parallel N / --resource_cap_pct N")
    print()

    # -- existing results --
    print(f"  EXISTING RESULTS (auto-detected, v5-compatible)")
    if cmd in ("all", "fastp") and params.get("use_fastp", True):
        print(f"    FASTP  : {fastp_done:>4d} / {total} completed"
              + (f"   <- will skip {fastp_done}" if fastp_done else ""))
    if cmd in ("all", "spades", "spades2busco"):
        print(f"    SPAdes : {spades_done:>4d} / {total} completed"
              + (f"   <- will skip {spades_done}" if spades_done else ""))
    if cmd in ("all", "busco", "spades2busco"):
        print(f"    BUSCO  : {busco_done:>4d} / {total} completed"
              + (f"   <- will skip {busco_done}" if busco_done else ""))
    print()

    # -- BUSCO fasta preview (first 5 samples) --
    if cmd in ("all", "busco", "spades2busco"):
        print(f"  BUSCO FASTA PREVIEW (first 5 samples)")
        for s in samples[:5]:
            fasta = choose_busco_fasta(s, out_root, params)
            exists = fasta.exists() and fasta.stat().st_size > 0
            status = "OK" if exists else "NOT FOUND"
            print(f"    {s.sid:<30s} -> {fasta.name}  [{status}]")
        if total > 5:
            print(f"    ... and {total - 5} more samples")
        print()

    # -- sample preview --
    print(f"  SAMPLE PREVIEW (first 5)")
    for s in samples[:5]:
        r1_short = Path(s.r1).name if s.r1 else "(none)"
        print(f"    {s.sid:<30s}  R1={r1_short}")
    if total > 5:
        print(f"    ... and {total - 5} more")
    print()

    print(f"{sep}")

    return {
        "fastp_done": fastp_done, "fastp_todo": fastp_todo,
        "spades_done": spades_done, "spades_todo": spades_todo,
        "busco_done": busco_done, "busco_todo": busco_todo,
    }


def ask_user_confirmation(auto_yes: bool) -> bool:
    """Ask user to confirm. Returns True to proceed, False to abort."""
    if auto_yes:
        print("  --yes flag detected, proceeding automatically.\n")
        return True

    print("  Proceed with this configuration?")
    print("    [Y] Yes, run pipeline")
    print("    [N] No, abort")
    print("    [D] Show sample details")
    print()
    try:
        while True:
            choice = input("  Your choice [Y/n/d]: ").strip().lower()
            if choice in ("", "y", "yes"):
                print()
                return True
            elif choice in ("n", "no", "q", "quit"):
                print("\n  Aborted by user. Adjust parameters and re-run.\n")
                return False
            elif choice in ("d", "detail", "details"):
                print("\n  (Details are shown above in the execution plan.)")
                print("  Re-run with adjusted parameters if needed.\n")
                continue
            else:
                print(f"  Unknown choice: '{choice}'. Enter Y, N, or D.")
    except (EOFError, KeyboardInterrupt):
        print("\n\n  Aborted by user.\n")
        return False


# ==============================================================
# 3. Sample dataclass
# ==============================================================

@dataclass(frozen=True)
class Sample:
    sid: str
    r1: str
    r2: str


# ==============================================================
# 4. R1/R2 suffix knowledge base + auto-detection (ENHANCED v5)
# ==============================================================

_SUFFIX_PAIRS: List[Tuple[str, str]] = [
    # Illumina / standard
    ("_R1_001.fastq.gz", "_R2_001.fastq.gz"),
    ("_R1_001.fq.gz", "_R2_001.fq.gz"),
    ("_R1_001.fastq", "_R2_001.fastq"),
    ("_R1_001.fq", "_R2_001.fq"),
    # Illumina lane patterns
    ("_L001_R1_001.fastq.gz", "_L001_R2_001.fastq.gz"),
    ("_L002_R1_001.fastq.gz", "_L002_R2_001.fastq.gz"),
    ("_L003_R1_001.fastq.gz", "_L003_R2_001.fastq.gz"),
    ("_L004_R1_001.fastq.gz", "_L004_R2_001.fastq.gz"),
    # Novogene / common Chinese sequencing providers
    (".R1.raw.fastq.gz", ".R2.raw.fastq.gz"),
    (".R1.clean.fastq.gz", ".R2.clean.fastq.gz"),
    ("_R1.raw.fastq.gz", "_R2.raw.fastq.gz"),
    ("_R1.raw.fq.gz", "_R2.raw.fq.gz"),
    ("_R1.clean.fq.gz", "_R2.clean.fq.gz"),
    # BGI / MGI
    ("_1.fq.gz", "_2.fq.gz"),
    ("_1.fastq.gz", "_2.fastq.gz"),
    ("_1.fq", "_2.fq"),
    ("_1.fastq", "_2.fastq"),
    ("_1.clean.fq.gz", "_2.clean.fq.gz"),
    ("_1.clean.fastq.gz", "_2.clean.fastq.gz"),
    # Generic _R1/_R2
    ("_R1.fastq.gz", "_R2.fastq.gz"),
    ("_R1.fq.gz", "_R2.fq.gz"),
    ("_R1.fastq", "_R2.fastq"),
    ("_R1.fq", "_R2.fq"),
    # Generic .R1/.R2
    (".R1.fastq.gz", ".R2.fastq.gz"),
    (".R1.fq.gz", ".R2.fq.gz"),
    (".R1.fastq", ".R2.fastq"),
    (".R1.fq", ".R2.fq"),
    # Macrogen
    ("_R1_trimmed.fastq.gz", "_R2_trimmed.fastq.gz"),
    ("_R1_trimmed.fq.gz", "_R2_trimmed.fq.gz"),
    # SRA / older
    ("_1.sanfastq.gz", "_2.sanfastq.gz"),
]

_SUFFIX_PRESETS: Dict[str, List[Tuple[str, str]]] = {
    "illumina": [
        ("_R1_001.fastq.gz", "_R2_001.fastq.gz"),
        ("_R1_001.fq.gz", "_R2_001.fq.gz"),
        ("_L001_R1_001.fastq.gz", "_L001_R2_001.fastq.gz"),
        ("_L002_R1_001.fastq.gz", "_L002_R2_001.fastq.gz"),
        ("_L003_R1_001.fastq.gz", "_L003_R2_001.fastq.gz"),
        ("_L004_R1_001.fastq.gz", "_L004_R2_001.fastq.gz"),
    ],
    "bgi": [
        ("_1.fq.gz", "_2.fq.gz"),
        ("_1.fastq.gz", "_2.fastq.gz"),
        ("_1.clean.fq.gz", "_2.clean.fq.gz"),
        ("_1.clean.fastq.gz", "_2.clean.fastq.gz"),
    ],
    "mgi": [
        ("_1.fq.gz", "_2.fq.gz"),
        ("_1.fastq.gz", "_2.fastq.gz"),
        ("_1.clean.fq.gz", "_2.clean.fq.gz"),
    ],
    "novogene": [
        (".R1.raw.fastq.gz", ".R2.raw.fastq.gz"),
        (".R1.clean.fastq.gz", ".R2.clean.fastq.gz"),
        ("_R1.raw.fastq.gz", "_R2.raw.fastq.gz"),
        ("_1.fq.gz", "_2.fq.gz"),
    ],
    "macrogen": [
        ("_R1_trimmed.fastq.gz", "_R2_trimmed.fastq.gz"),
        ("_R1_001.fastq.gz", "_R2_001.fastq.gz"),
    ],
    "generic": [
        ("_R1.fastq.gz", "_R2.fastq.gz"),
        ("_R1.fq.gz", "_R2.fq.gz"),
        (".R1.fastq.gz", ".R2.fastq.gz"),
        (".R1.fq.gz", ".R2.fq.gz"),
        ("_1.fastq.gz", "_2.fastq.gz"),
        ("_1.fq.gz", "_2.fq.gz"),
    ],
}


def resolve_suffix_fq(raw_args: List[str]) -> List[Tuple[str, str]]:
    if not raw_args:
        return list(_SUFFIX_PAIRS)
    if len(raw_args) == 1 and raw_args[0].lower() == "auto":
        return list(_SUFFIX_PAIRS)
    all_presets = all(a.lower() in _SUFFIX_PRESETS for a in raw_args)
    if all_presets:
        pairs: List[Tuple[str, str]] = []
        seen: Set[Tuple[str, str]] = set()
        for name in raw_args:
            for pair in _SUFFIX_PRESETS[name.lower()]:
                if pair not in seen:
                    seen.add(pair)
                    pairs.append(pair)
        return pairs
    preset_names = [a for a in raw_args if a.lower() in _SUFFIX_PRESETS]
    non_preset = [a for a in raw_args if a.lower() not in _SUFFIX_PRESETS
                  and a.lower() != "auto"]
    if preset_names and non_preset and len(non_preset) % 2 == 0:
        pairs = []
        seen_set: Set[Tuple[str, str]] = set()
        for name in preset_names:
            for pair in _SUFFIX_PRESETS[name.lower()]:
                if pair not in seen_set:
                    seen_set.add(pair)
                    pairs.append(pair)
        for i in range(0, len(non_preset), 2):
            pair = (non_preset[i], non_preset[i + 1])
            if pair not in seen_set:
                seen_set.add(pair)
                pairs.append(pair)
        return pairs
    if len(non_preset) % 2 == 0 and len(non_preset) >= 2 and not preset_names:
        pairs = []
        for i in range(0, len(non_preset), 2):
            pairs.append((non_preset[i], non_preset[i + 1]))
        return pairs
    raise ValueError(
        f"Invalid --suffix_fq: {' '.join(raw_args)}\n"
        f"\n"
        f"  Expected one of:\n"
        f"    --suffix_fq auto                            (detect automatically)\n"
        f"    --suffix_fq illumina                        (vendor preset)\n"
        f"    --suffix_fq illumina bgi                    (combine presets)\n"
        f"    --suffix_fq _R1.fastq.gz _R2.fastq.gz      (explicit pair)\n"
        f"    --suffix_fq _R1.fq.gz _R2.fq.gz  _1.fq.gz _2.fq.gz  (multiple pairs)\n"
        f"\n"
        f"  Available presets: {', '.join(sorted(_SUFFIX_PRESETS.keys()))}\n"
        f"\n"
        f"  Note: explicit suffixes must come in R1/R2 pairs (even number of args)")


def _swap_r1_to_r2(r1_path: str, r1_suf: str, r2_suf: str) -> Optional[str]:
    if r1_path.endswith(r1_suf):
        candidate = r1_path[: -len(r1_suf)] + r2_suf
        if Path(candidate).exists():
            return candidate
    return None


def auto_find_r2(r1: str,
                 pairs: Optional[List[Tuple[str, str]]] = None) -> Optional[str]:
    for r1s, r2s in (pairs or _SUFFIX_PAIRS):
        r2 = _swap_r1_to_r2(r1, r1s, r2s)
        if r2:
            return r2
    return None


def strip_r1_suffix(filename: str,
                    pairs: Optional[List[Tuple[str, str]]] = None) -> str:
    for r1s, _ in (pairs or _SUFFIX_PAIRS):
        if filename.endswith(r1s):
            return filename[: -len(r1s)]
    return re.sub(r"\.(fastq|fq)(\.gz)?$", "", filename, flags=re.I)


def auto_detect_suffix_in_dir(
    dir_path: Path,
    pairs: Optional[List[Tuple[str, str]]] = None,
) -> List[Tuple[str, str]]:
    if not dir_path.is_dir():
        return []
    fnames = {f.name for f in dir_path.iterdir() if f.is_file()}
    if not fnames:
        return []
    hits: List[Tuple[int, Tuple[str, str]]] = []
    for r1s, r2s in (pairs or _SUFFIX_PAIRS):
        r1_files = [f for f in fnames if f.endswith(r1s) and len(f) > len(r1s)]
        if not r1_files:
            continue
        n_paired = sum(1 for f in r1_files if (f[: -len(r1s)] + r2s) in fnames)
        if n_paired > 0:
            hits.append((n_paired, (r1s, r2s)))
    hits.sort(key=lambda x: -x[0])
    return [pair for _, pair in hits]


def detect_pairs_in_dir(dir_path: Path, r1_suf: str,
                        r2_suf: str) -> Dict[str, Tuple[Path, Path]]:
    pairs: Dict[str, Tuple[Path, Path]] = {}
    if not dir_path.is_dir():
        return pairs
    for fq1 in sorted(dir_path.iterdir()):
        if not fq1.is_file():
            continue
        name = fq1.name
        if not name.endswith(r1_suf) or len(name) <= len(r1_suf):
            continue
        sample = name[: -len(r1_suf)]
        fq2 = dir_path / f"{sample}{r2_suf}"
        if fq2.is_file():
            pairs[sample] = (fq1.resolve(), fq2.resolve())
    return pairs


def detect_pairs_multi_suffix(
    dir_path: Path,
    suffix_pairs: List[Tuple[str, str]],
) -> Dict[str, Tuple[Path, Path]]:
    result: Dict[str, Tuple[Path, Path]] = {}
    for r1s, r2s in suffix_pairs:
        found = detect_pairs_in_dir(dir_path, r1s, r2s)
        for sid, paths in found.items():
            if sid not in result:
                result[sid] = paths
    return result


# ==============================================================
# 5. Sample table parser
# ==============================================================

_COL_SAMPLE = {"sample", "name", "sample_name", "id", "sample_id"}
_COL_SDIR = {"sample_dir", "subdir", "dir", "directory", "outdir"}
_COL_FQ1 = {"fq1", "r1", "fastq1", "abs_fastq", "abs_fq1",
             "read1", "reads1", "fq_r1", "fq_1"}
_COL_FQ2 = {"fq2", "r2", "fastq2", "abs_fq2",
             "read2", "reads2", "fq_r2", "fq_2"}


def _norm(s: str) -> str:
    return s.lstrip("#").strip().lower()


def _looks_like_fastq(s: str) -> bool:
    return bool(re.search(r"\.(fastq|fq)(\.gz)?$", s, re.I))


def _looks_like_abs_dir(s: str) -> bool:
    p = Path(s)
    return p.is_absolute() and p.is_dir()


def parse_sample_table(
    table_path: Path,
    suffix_pairs: Optional[List[Tuple[str, str]]] = None,
    require_fastq: bool = True,
) -> Tuple[List[Sample], List[str]]:
    lines = table_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    spairs = suffix_pairs or _SUFFIX_PAIRS
    header_idx = -1
    col_map: Dict[str, int] = {}
    for i, raw_line in enumerate(lines):
        stripped = raw_line.strip()
        if not stripped:
            continue
        # FIX v6.2: also detect non-commented header lines (tab/space-delimited
        # files where the header has no leading '#').  Previously only '#'-prefixed
        # headers were recognised; bare headers were parsed as fake sample rows,
        # emitting Sample(sid='sample', r1='r1', r2='r2').
        candidate = stripped.lstrip("#")
        parts = re.split(r"[\t ]+", candidate)
        norms = [_norm(p) for p in parts]
        all_known = (_COL_SAMPLE | _COL_SDIR | _COL_FQ1 | _COL_FQ2)
        if any(n in all_known for n in norms):
            header_idx = i
            for ci, n in enumerate(norms):
                if n in _COL_SAMPLE and "sample" not in col_map:
                    col_map["sample"] = ci
                elif n in _COL_SDIR and "sdir" not in col_map:
                    col_map["sdir"] = ci
                elif n in _COL_FQ1 and "fq1" not in col_map:
                    col_map["fq1"] = ci
                elif n in _COL_FQ2 and "fq2" not in col_map:
                    col_map["fq2"] = ci
            break
        if stripped.startswith("#"):
            # Commented line that didn't match any known header — skip
            continue
    samples: List[Sample] = []
    warnings: List[str] = []
    seen_sids: Set[str] = set()
    data_lines = lines[header_idx + 1:] if header_idx >= 0 else lines
    for line_offset, raw_line in enumerate(data_lines):
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        lineno = header_idx + 2 + line_offset
        parts = re.split(r"\t", stripped)
        if len(parts) == 1:
            parts = re.split(r" {2,}", stripped)
        if len(parts) == 1:
            parts = stripped.split()
        sid: Optional[str] = None
        sdir: Optional[str] = None
        fq1: Optional[str] = None
        fq2: Optional[str] = None
        if col_map:
            def _get(role: str) -> Optional[str]:
                ci = col_map.get(role)
                if ci is not None and ci < len(parts):
                    v = parts[ci].strip()
                    return v if v else None
                return None
            sid = _get("sample")
            sdir = _get("sdir")
            fq1 = _get("fq1")
            fq2 = _get("fq2")
        else:
            fq_cols = [i for i, p in enumerate(parts) if _looks_like_fastq(p)]
            dir_cols = [i for i, p in enumerate(parts)
                        if _looks_like_abs_dir(p) and not _looks_like_fastq(p)]
            other_cols = [i for i in range(len(parts))
                          if i not in fq_cols and i not in dir_cols]
            if len(fq_cols) >= 2:
                fq1, fq2 = parts[fq_cols[0]], parts[fq_cols[1]]
            elif len(fq_cols) == 1:
                fq1 = parts[fq_cols[0]]
            if dir_cols:
                sdir = parts[dir_cols[0]]
            if other_cols:
                sid = parts[other_cols[0]]
        if sdir and not fq1:
            sdir_path = Path(sdir)
            if sdir_path.is_dir():
                detected = auto_detect_suffix_in_dir(sdir_path, spairs)
                if detected:
                    best_r1s, best_r2s = detected[0]
                    cands = sorted(sdir_path.glob(f"*{best_r1s}"))
                    if cands:
                        fq1 = str(cands[0])
                        r2c = sdir_path / (cands[0].name[: -len(best_r1s)] + best_r2s)
                        if r2c.exists():
                            fq2 = str(r2c)
        if fq1 and not fq2:
            fq2 = auto_find_r2(fq1, spairs)
        if not sid and fq1:
            sid = strip_r1_suffix(Path(fq1).name, spairs)
        if not sid and not require_fastq:
            pass
        elif not sid:
            sid = f"sample_line{lineno}"
        if not sid:
            warnings.append(f"Line {lineno}: cannot determine sample name -- skipped")
            continue
        sid = safe_name(sid)
        if require_fastq:
            if not fq1:
                warnings.append(f"Line {lineno} [{sid}]: R1 not found -- skipped")
                continue
            if not fq2:
                warnings.append(
                    f"Line {lineno} [{sid}]: R2 not found "
                    f"(R1={fq1}; tried {len(spairs)} patterns) -- skipped")
                continue
            fq1_p, fq2_p = Path(fq1), Path(fq2)
            if not fq1_p.exists():
                warnings.append(
                    f"Line {lineno} [{sid}]: R1 does not exist: {fq1} -- skipped")
                continue
            if not fq2_p.exists():
                warnings.append(
                    f"Line {lineno} [{sid}]: R2 does not exist: {fq2} -- skipped")
                continue
            fq1 = str(fq1_p.resolve())
            fq2 = str(fq2_p.resolve())
        else:
            fq1 = str(Path(fq1).resolve()) if fq1 and Path(fq1).exists() else ""
            fq2 = str(Path(fq2).resolve()) if fq2 and Path(fq2).exists() else ""
        base_sid, cnt = sid, 2
        while sid in seen_sids:
            sid = f"{base_sid}_{cnt}"
            cnt += 1
        if sid != base_sid:
            warnings.append(f"Line {lineno}: duplicate '{base_sid}' -> '{sid}'")
        seen_sids.add(sid)
        samples.append(Sample(sid=sid, r1=fq1 or "", r2=fq2 or ""))
    return samples, warnings


# ==============================================================
# 6. Directory-based sample scanning
# ==============================================================

def _collect_tasks_from_dir(
    input_dir: Path,
    fq_position: str,
    subdir_list: Optional[List[str]],
    suffix_pairs: List[Tuple[str, str]],
    exclude_dirs: Optional[Set[Path]] = None,
) -> Tuple[List[Tuple[str, Optional[Path], Optional[Path]]],
           List[Tuple[str, str]]]:
    tasks: List[Tuple[str, Optional[Path], Optional[Path]]] = []
    used_suffixes: List[Tuple[str, str]] = []
    _excl = exclude_dirs or set()
    if fq_position == "flat":
        merged = detect_pairs_multi_suffix(input_dir, suffix_pairs)
        allow = set(subdir_list) if subdir_list else None
        for sample in sorted(merged):
            if allow and sample not in allow:
                continue
            fq1, fq2 = merged[sample]
            tasks.append((sample, fq1, fq2))
        detected = auto_detect_suffix_in_dir(input_dir, suffix_pairs)
        used_suffixes = detected[:5]
        return tasks, used_suffixes
    groups = subdir_list or sorted(
        p.name for p in input_dir.iterdir()
        if p.is_dir() and p.resolve() not in _excl)
    for g in groups:
        gdir = input_dir / g
        if not gdir.is_dir():
            logging.warning(f"[SCAN] subdir not found: {gdir}")
            continue
        if gdir.resolve() in _excl:
            logging.debug(f"[SCAN] skipping output dir: {gdir}")
            continue
        merged = detect_pairs_multi_suffix(gdir, suffix_pairs)
        if not merged:
            tasks.append((g, None, None))
        elif len(merged) == 1:
            _, (fq1, fq2) = next(iter(merged.items()))
            tasks.append((g, fq1, fq2))
        else:
            for samp, (fq1, fq2) in sorted(merged.items()):
                tasks.append((f"{g}__{samp}", fq1, fq2))
        detected = auto_detect_suffix_in_dir(gdir, suffix_pairs)
        for d in detected:
            if d not in used_suffixes:
                used_suffixes.append(d)
    return tasks, used_suffixes[:5]


def collect_samples_from_dir(
    args,
    suffix_pairs: List[Tuple[str, str]],
) -> Tuple[List[Sample], List[str]]:
    input_dir = Path(args.input_dir).resolve()
    if not input_dir.is_dir():
        # FIX v6.2: raise instead of sys.exit() — sys.exit() inside
        # a ThreadPoolExecutor worker only kills that thread, not the
        # process, causing silent data loss.  Let callers handle it.
        raise RuntimeError(f"[INPUT] --input_dir not found: {input_dir}")
    subdir_list = resolve_list_arg(
        getattr(args, "subdir_list", None), input_dir)
    # v6.16: skip the output directory if it lives inside the input directory
    excl: Set[Path] = set()
    out_raw = getattr(args, "out_root", None)
    if out_raw:
        out_p = Path(out_raw).resolve()
        # Add the output root and its immediate children (fastp_dir, etc.)
        excl.add(out_p)
    tasks, used_suffixes = _collect_tasks_from_dir(
        input_dir, args.fq_position, subdir_list, suffix_pairs,
        exclude_dirs=excl)
    if used_suffixes:
        suf_str = "  ".join(f"{r1} / {r2}" for r1, r2 in used_suffixes[:3])
        logging.info(f"[SCAN] detected suffix pattern(s): {suf_str}")
    logging.info(
        f"[SCAN] input_dir={input_dir}  mode={args.fq_position}  "
        f"tasks={len(tasks)}  suffix_patterns_tried={len(suffix_pairs)}")
    samples: List[Sample] = []
    missing: List[str] = []
    for tid, fq1, fq2 in tasks:
        sid = safe_name(tid)
        if fq1 is None or fq2 is None:
            missing.append(sid)
        else:
            samples.append(Sample(sid=sid, r1=str(fq1), r2=str(fq2)))
    return samples, sorted(set(missing))


# ==============================================================
# 7. Unified sample collection
# ==============================================================

def collect_samples(
    args,
    suffix_pairs: List[Tuple[str, str]],
    require_fastq: bool = True,
) -> Tuple[List[Sample], List[str]]:
    if getattr(args, "sample_table", None):
        tbl = Path(args.sample_table)
        if not tbl.is_file():
            # FIX v6.2: raise instead of sys.exit() for same reason as above
            raise RuntimeError(f"[INPUT] sample_table not found: {tbl}")
        samples, warns = parse_sample_table(
            tbl, suffix_pairs=suffix_pairs, require_fastq=require_fastq)
        for w in warns:
            logging.warning(f"[TABLE] {w}")
        logging.info(f"[TABLE] {tbl}  parsed={len(samples)} samples")
        return samples, []
    if getattr(args, "input_dir", None):
        return collect_samples_from_dir(args, suffix_pairs)
    return [], []


# ==============================================================
# 8. Comprehensive pre-flight validation (v6.3)
#
#   Design principle: collect ALL errors/warnings in one pass,
#   report everything at once, then HARD-STOP if any error exists.
#   Nothing runs until every check is GREEN.
#
#   Check categories (in order):
#     A. Python version
#     B. Conda / mamba availability
#     C. Tool binaries (fastp / spades / busco) — version probe
#     D. Snakemake (if --smk)
#     E. System resources (CPUs, RAM, disk)
#     F. Parameter sanity (threads, memory, modes, flags)
#     G. Sample integrity — ALL samples, no cap:
#          G1. R1/R2 file existence
#          G2. R1/R2 file size > 0
#          G3. R1/R2 gzip magic bytes (first 2 bytes = 0x1f 0x8b)
#          G4. R1/R2 are not the same file (symlink / identical path)
#          G5. Sample ID valid (no shell-unsafe chars)
#          G6. No duplicate sample IDs
#     H. Output directory writable
#     I. BUSCO lineage path (if local) or reachability hint
#     J. BUSCO fasta_dir (if --busco_fasta_dir)
#     K. For spades2busco / all: fastp output pre-check (if use_fastp)
#     L. SPAdes memory <= system RAM
#     M. max_parallel * threads_spades <= system CPUs (warning only)
# ==============================================================

# ── internal result carrier ──────────────────────────────────────

class _Check:
    """A single preflight result."""
    __slots__ = ("category", "label", "ok", "message", "fatal")

    def __init__(self, category: str, label: str,
                 ok: bool, message: str, fatal: bool = True):
        self.category = category
        self.label    = label
        self.ok       = ok
        self.message  = message
        self.fatal    = fatal   # False → warning only, does not block run


def _check_gzip_magic(path: str) -> bool:
    """Return True if *path* starts with gzip magic bytes 0x1f 0x8b."""
    try:
        with open(path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except Exception:
        return False


def _is_plain_fastq(path: str) -> bool:
    """Return True if path looks like uncompressed fastq (starts with '@')."""
    try:
        with open(path, "rb") as f:
            first = f.read(1)
            return first == b"@"
    except Exception:
        return False


# ── A. Python version ────────────────────────────────────────────

def _chk_python() -> _Check:
    pv = sys.version_info
    ok = pv >= (3, 7)
    ver = f"{pv.major}.{pv.minor}.{pv.micro}"
    msg = (f"Python {ver}  ({sys.executable})"
           if ok else
           f"Python {ver} is too old — need >= 3.7  ({sys.executable})")
    return _Check("A", "Python version", ok, msg)


# ── B. Conda availability ────────────────────────────────────────

def _chk_conda(conda_exe: str, needs_conda: bool) -> Optional[_Check]:
    if not needs_conda:
        return None
    found = which_or_none(conda_exe)
    if found:
        return _Check("B", f"conda ({conda_exe})", True,
                      f"found at {found}")
    return _Check("B", f"conda ({conda_exe})", False,
                  f"'{conda_exe}' not found in PATH — required because "
                  f"--conda_env_* is set.\n"
                  f"  Fix: conda activate <env>  or  export PATH=...:$PATH")


# ── C. Tool binary probes ────────────────────────────────────────

def _probe_tool(label: str, exe_name: str,
                conda_env: Optional[str], tool_path: Optional[str],
                conda_exe: str) -> _Check:
    cat = "C"
    try:
        cmd = build_tool_cmd(
            tool_label=label, exe_name=exe_name,
            base_args=["--version"],
            conda_env=conda_env, tool_path=tool_path,
            conda_exe=conda_exe)
    except FileNotFoundError as e:
        return _Check(cat, label, False,
                      f"binary not found — {e}\n"
                      f"  Fix: conda activate <env>  or  "
                      f"--{label.lower()}_path /path/to/{exe_name}")
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, text=True, timeout=30)
        out_txt = (p.stdout + p.stderr).strip()
        ver_line = out_txt.splitlines()[0] if out_txt else "(no output)"
        ok = p.returncode == 0
        source_bits = []
        if conda_env:
            source_bits.append(f"env={conda_env}")
        if tool_path:
            source_bits.append(f"path={tool_path}")
        source_txt = ("  " + "; ".join(source_bits)) if source_bits else ""
        msg = (f"cmd={' '.join(shlex.quote(str(x)) for x in cmd)}{source_txt}  version: {ver_line}"
               if ok else
               f"cmd={' '.join(shlex.quote(str(x)) for x in cmd)}{source_txt}  --version exited rc={p.returncode}: {ver_line}")
        return _Check(cat, label, ok, msg, fatal=not ok)
    except subprocess.TimeoutExpired:
        return _Check(cat, label, True,
                      f"cmd={' '.join(shlex.quote(str(x)) for x in cmd)}  (--version timed out — assumed OK)")
    except Exception as e:
        return _Check(cat, label, False, f"error probing tool: {e}")


# ── D. Snakemake ─────────────────────────────────────────────────

def _chk_snakemake(use_smk: bool) -> Optional[_Check]:
    if not use_smk:
        return None
    found = which_or_none("snakemake")
    if found:
        return _Check("D", "snakemake", True, f"found at {found}")
    return _Check("D", "snakemake", False,
                  "snakemake not found in PATH but --smk was requested.\n"
                  "  Fix: pip install snakemake  or  "
                  "conda install -c bioconda snakemake")


# ── E. System resources ──────────────────────────────────────────

def _chk_resources(params: dict, out_root: Path) -> List[_Check]:
    results: List[_Check] = []
    cpus, mem_gb = get_system_resources()
    results.append(_Check("E", "System CPUs/RAM", True,
                          f"CPUs={cpus}  RAM={mem_gb} GB"))

    # Disk space
    free_gb = get_disk_free_gb(out_root)
    if free_gb >= 0:
        disk_ok = free_gb >= 5.0
        results.append(_Check(
            "E", "Disk free space",
            disk_ok,
            f"{free_gb:.1f} GB free on {out_root.anchor}"
            + ("" if disk_ok else
               " — WARNING: < 5 GB; assembly may run out of space"),
            fatal=False))   # warning only; some pipelines need very little space

    # SPAdes memory vs system RAM
    spades_mem = int(params.get("spades_mem_gb", 0) or 0)
    if spades_mem > 0 and spades_mem > mem_gb:
        results.append(_Check(
            "E", "SPAdes memory vs RAM", False,
            f"--spades_mem_gb {spades_mem} exceeds system RAM {mem_gb} GB.\n"
            f"  Fix: --spades_mem_gb {max(1, mem_gb - 2)}"))

    # max_parallel * threads_spades vs CPUs (warning)
    par = int(params.get("max_parallel", 8) or 8)
    thr = int(params.get("threads_spades", 16) or 16)
    if par * thr > cpus * 2:
        results.append(_Check(
            "E", "Parallelism vs CPUs", False,
            f"max_parallel({par}) × threads_spades({thr}) = {par*thr} "
            f"greatly exceeds CPUs ({cpus}).\n"
            f"  Suggestion: --max_parallel {max(1, cpus // thr)} "
            f"or --threads_spades {max(1, cpus // par)}",
            fatal=False))
    return results


# ── F. Parameter sanity ──────────────────────────────────────────

def _chk_params(cmd: str, params: dict) -> List[_Check]:
    results: List[_Check] = []

    # Thread counts
    for key in ("threads_fastp", "threads_spades", "threads_busco"):
        val = int(params.get(key, 0) or 0)
        if val < 1:
            results.append(_Check("F", key, False,
                                  f"--{key} = {val}; must be >= 1"))

    # Memory
    mem = int(params.get("spades_mem_gb", 0) or 0)
    if mem < 1 and cmd in ("all", "spades", "spades2busco"):
        results.append(_Check("F", "spades_mem_gb", False,
                               f"--spades_mem_gb = {mem}; must be >= 1"))

    # Assembly mode
    asm_mode = normalize_assembly_mode(params.get("assembly_mode", ""))
    if asm_mode and asm_mode not in VALID_ASM_MODES:
        results.append(_Check(
            "F", "assembly_mode", False,
            f"--assembly_mode '{asm_mode}' is invalid.\n"
            f"  Valid: {', '.join(sorted(VALID_ASM_MODES))}"))

    # BUSCO lineage required
    if cmd in ("all", "busco", "spades2busco"):
        if not params.get("busco_lineage"):
            results.append(_Check(
                "F", "busco_lineage", False,
                "--busco_lineage is required.\n"
                "  Example: --busco_lineage metazoa_odb10\n"
                "  List available: busco --list-datasets"))
        else:
            # Lineage is local path but does not exist
            lin = params["busco_lineage"]
            if is_local_busco_lineage(lin) and not Path(lin).exists():
                results.append(_Check(
                    "F", "busco_lineage path", False,
                    f"--busco_lineage looks like a local path but was not found:\n"
                    f"  {lin}"))

    # suffix_spades_out shell-safety
    suffix = params.get("suffix_spades_out", "")
    if suffix and re.search(r'[/\\:*?"<>| ]', suffix):
        results.append(_Check(
            "F", "suffix_spades_out", False,
            f"--suffix_spades_out '{suffix}' contains shell-unsafe characters.\n"
            f"  Use only: letters, digits, underscore, hyphen, dot."))

    # resource_cap_pct range
    cap = float(params.get("resource_cap_pct", 90) or 90)
    if not (10 <= cap <= 100):
        results.append(_Check(
            "F", "resource_cap_pct", False,
            f"--resource_cap_pct {cap} out of range 10–100."))

    # max_parallel >= 1
    par = int(params.get("max_parallel", 1) or 1)
    if par < 1:
        results.append(_Check("F", "max_parallel", False,
                               f"--max_parallel {par}; must be >= 1"))

    return results


# ── G. Sample integrity — ALL samples, zero cap ──────────────────

def _chk_samples(cmd: str, samples: List["Sample"],
                 params: dict) -> List[_Check]:
    """Check every sample for path existence, size, gzip magic, duplicates."""
    results: List[_Check] = []
    if not samples:
        results.append(_Check("G", "sample count", False,
                               "No samples found — check -i / --sample_table."))
        return results

    # G6: duplicate sample IDs
    seen_ids: Dict[str, int] = {}
    for s in samples:
        seen_ids[s.sid] = seen_ids.get(s.sid, 0) + 1
    dups = [sid for sid, cnt in seen_ids.items() if cnt > 1]
    if dups:
        results.append(_Check(
            "G", "duplicate sample IDs", False,
            f"{len(dups)} duplicate sample ID(s): "
            f"{', '.join(dups[:10])}"
            + ("..." if len(dups) > 10 else "")))

    # G5: sample ID shell safety
    bad_ids = [s.sid for s in samples
               if re.search(r"[^\w.\-+]", s.sid)]
    if bad_ids:
        results.append(_Check(
            "G", "sample ID characters", False,
            f"{len(bad_ids)} sample ID(s) contain shell-unsafe characters "
            f"(spaces, quotes, $, etc.):\n"
            f"  {', '.join(bad_ids[:10])}"
            + ("..." if len(bad_ids) > 10 else "")))

    # G1–G4: per-file checks (only for commands that need FASTQ)
    if cmd not in ("all", "fastp", "spades", "spades2busco"):
        return results

    missing_r1:  List[str] = []
    missing_r2:  List[str] = []
    empty_r1:    List[str] = []
    empty_r2:    List[str] = []
    notgz_r1:    List[str] = []
    notgz_r2:    List[str] = []
    same_file:   List[str] = []

    for s in samples:
        if not s.r1 and not s.r2:
            continue  # BUSCO-only samples without fastq

        r1p = Path(s.r1) if s.r1 else None
        r2p = Path(s.r2) if s.r2 else None

        # G1: existence
        if r1p is None or not r1p.exists():
            missing_r1.append(f"  [{s.sid}] {s.r1}")
            continue   # skip further checks for this sample
        if r2p is None or not r2p.exists():
            missing_r2.append(f"  [{s.sid}] {s.r2}")
            continue

        # G2: size > 0
        if r1p.stat().st_size == 0:
            empty_r1.append(f"  [{s.sid}] {s.r1}")
        if r2p.stat().st_size == 0:
            empty_r2.append(f"  [{s.sid}] {s.r2}")

        # G3: gzip magic (only if .gz extension)
        if s.r1.endswith(".gz") and not _check_gzip_magic(s.r1):
            notgz_r1.append(f"  [{s.sid}] {s.r1}  (not a valid gzip file)")
        if s.r2.endswith(".gz") and not _check_gzip_magic(s.r2):
            notgz_r2.append(f"  [{s.sid}] {s.r2}  (not a valid gzip file)")

        # G4: R1 and R2 must not be the same file
        try:
            if r1p.resolve() == r2p.resolve():
                same_file.append(f"  [{s.sid}] R1 == R2: {s.r1}")
        except Exception:
            pass

    def _add(label: str, lst: List[str], hint: str = "") -> None:
        if lst:
            body = "\n".join(lst[:15])
            if len(lst) > 15:
                body += f"\n  ... and {len(lst) - 15} more"
            results.append(_Check(
                "G", label, False,
                f"{len(lst)} file(s) with issue:\n{body}"
                + (f"\n  Hint: {hint}" if hint else "")))

    _add("R1 files not found", missing_r1,
         "check --input_dir / --sample_table path and --suffix_fq pattern")
    _add("R2 files not found", missing_r2,
         "auto-detect failed; specify --suffix_fq explicitly")
    _add("R1 files are empty (0 bytes)", empty_r1)
    _add("R2 files are empty (0 bytes)", empty_r2)
    _add("R1 files have bad gzip header", notgz_r1,
         "file may be truncated or not actually gzip-compressed")
    _add("R2 files have bad gzip header", notgz_r2,
         "file may be truncated or not actually gzip-compressed")
    _add("R1 and R2 are the same file", same_file,
         "check your sample table for copy-paste errors")

    return results


# ── H. Output directory ──────────────────────────────────────────

def _chk_outdir(out_root: Path) -> _Check:
    try:
        out_root.mkdir(parents=True, exist_ok=True)
        test_f = out_root / ".preflight_write_test"
        test_f.write_text("ok", encoding="utf-8")
        test_f.unlink()
        return _Check("H", "output dir writable", True, str(out_root))
    except Exception as e:
        return _Check("H", "output dir writable", False,
                      f"Cannot write to output dir: {out_root}\n  {e}")


# ── I. BUSCO fasta_dir ───────────────────────────────────────────

def _chk_busco_fasta_dir(cmd: str, params: dict) -> Optional[_Check]:
    if cmd not in ("all", "busco", "spades2busco"):
        return None
    bfd = params.get("busco_fasta_dir")
    if not bfd:
        return None
    p = Path(bfd)
    if not p.is_dir():
        return _Check("I", "busco_fasta_dir", False,
                      f"--busco_fasta_dir not found: {bfd}")
    fa_files = (list(p.glob("*.fasta")) + list(p.glob("*.fa"))
                + list(p.glob("*.fna")))
    if not fa_files:
        return _Check("I", "busco_fasta_dir contents", False,
                      f"--busco_fasta_dir '{bfd}' contains no .fasta/.fa/.fna files.")
    return _Check("I", "busco_fasta_dir", True,
                  f"{bfd}  ({len(fa_files)} fasta files found)")


# ── J. Snakemake backend pre-check ───────────────────────────────

def _chk_smk_conda_envs(use_smk: bool, params: dict) -> List[_Check]:
    """When using snakemake + conda envs, verify the envs actually exist."""
    if not use_smk:
        return []
    results: List[_Check] = []
    conda_exe = params.get("conda_exe", "conda")
    for key, label in [("conda_env_fastp", "conda_env_fastp"),
                       ("conda_env_spades", "conda_env_spades"),
                       ("conda_env_busco", "conda_env_busco")]:
        env = params.get(key)
        if not env:
            continue
        try:
            p = subprocess.run(
                [conda_exe, "env", "list"],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, timeout=20)
            env_names = [ln.split()[0] for ln in p.stdout.splitlines()
                         if ln and not ln.startswith("#") and ln.split()]
            if env in env_names:
                results.append(_Check("J", f"conda env '{env}'", True,
                                      f"exists ({label})"))
            else:
                results.append(_Check(
                    "J", f"conda env '{env}'", False,
                    f"conda environment '{env}' not found "
                    f"(used by --{label}).\n"
                    f"  Available: {', '.join(env_names[:10])}\n"
                    f"  Fix: conda create -n {env} ..."))
        except Exception as e:
            results.append(_Check(
                "J", f"conda env '{env}'", False,
                f"Could not list conda envs: {e}"))
    return results


def _chk_smk_tool_env_binding(use_smk: bool, params: dict) -> List[_Check]:
    """Under Snakemake, verify each requested conda env can actually resolve its tool."""
    if not use_smk:
        return []
    results: List[_Check] = []
    tool_specs = [
        ("fastp", params.get("fastp_exe", "fastp"), params.get("conda_env_fastp"), params.get("fastp_path")),
        ("spades", params.get("spades_exe", "spades.py"), params.get("conda_env_spades"), params.get("spades_path")),
        ("busco", params.get("busco_exe", "busco"), params.get("conda_env_busco"), params.get("busco_path")),
    ]
    conda_exe = params.get("conda_exe", "conda")
    for label, exe_name, conda_env, tool_path in tool_specs:
        if not conda_env:
            continue
        try:
            cmd = build_tool_cmd(
                tool_label=label, exe_name=exe_name, base_args=["--help"],
                conda_env=conda_env, tool_path=tool_path, conda_exe=conda_exe)
            p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, timeout=40)
            combined = (p.stdout + p.stderr).strip()
            line = combined.splitlines()[0] if combined else "(no output)"
            ok = p.returncode in (0, 1)
            results.append(_Check(
                "J", f"{label} env binding", ok,
                (f"env={conda_env} resolves {exe_name}: {line}" if ok else
                 f"env={conda_env} could not run {exe_name} (rc={p.returncode}): {line}"),
                fatal=not ok))
        except Exception as e:
            results.append(_Check("J", f"{label} env binding", False,
                                  f"env={conda_env} probe failed for {exe_name}: {e}", fatal=True))
    return results


# ── Master preflight function ────────────────────────────────────

def preflight_full_check(cmd: str, params: dict,
                         samples: List["Sample"],
                         out_root: Path,
                         use_smk: bool = False) -> int:
    """
    Run ALL pre-flight checks.  Print a structured report.
    Return 0 only when every FATAL check passes.
    Warnings (fatal=False) are printed but do not block execution.
    """
    all_checks: List[_Check] = []
    conda_exe = params.get("conda_exe", "conda")

    # A. Python
    all_checks.append(_chk_python())

    # B. Conda
    needs_conda = any(params.get(k) for k in
                      ("conda_env_fastp", "conda_env_spades", "conda_env_busco"))
    chk_b = _chk_conda(conda_exe, needs_conda)
    if chk_b:
        all_checks.append(chk_b)

    # C. Tool binaries
    if cmd in ("all", "fastp") and params.get("use_fastp", True):
        all_checks.append(_probe_tool(
            "fastp", params.get("fastp_exe", "fastp"),
            params.get("conda_env_fastp"), params.get("fastp_path"), conda_exe))
    if cmd in ("all", "spades", "spades2busco"):
        all_checks.append(_probe_tool(
            "spades", params.get("spades_exe", "spades.py"),
            params.get("conda_env_spades"), params.get("spades_path"), conda_exe))
    if cmd in ("all", "busco", "spades2busco"):
        all_checks.append(_probe_tool(
            "busco", params.get("busco_exe", "busco"),
            params.get("conda_env_busco"), params.get("busco_path"), conda_exe))

    # D. Snakemake
    chk_d = _chk_snakemake(use_smk)
    if chk_d:
        all_checks.append(chk_d)

    # E. System resources
    all_checks.extend(_chk_resources(params, out_root))

    # F. Parameter sanity
    all_checks.extend(_chk_params(cmd, params))

    # G. Sample integrity (all samples, no cap)
    all_checks.extend(_chk_samples(cmd, samples, params))

    # H. Output dir
    all_checks.append(_chk_outdir(out_root))

    # I. BUSCO fasta dir
    chk_i = _chk_busco_fasta_dir(cmd, params)
    if chk_i:
        all_checks.append(chk_i)

    # J. Conda env existence (snakemake mode)
    all_checks.extend(_chk_smk_conda_envs(use_smk, params))
    all_checks.extend(_chk_smk_tool_env_binding(use_smk, params))

    # ── Print report ─────────────────────────────────────────────
    sep = "=" * 66
    print(f"\n{sep}")
    print(f"  PRE-FLIGHT CHECK REPORT  (v6.25 — {len(all_checks)} checks)")
    print(f"{sep}")

    errors:   List[_Check] = []
    warnings: List[_Check] = []
    passed:   List[_Check] = []

    for chk in all_checks:
        if chk.ok:
            passed.append(chk)
        elif chk.fatal:
            errors.append(chk)
        else:
            warnings.append(chk)

    # Always show errors and warnings; show passed only in verbose
    if errors:
        print(f"\n  ✗ ERRORS  ({len(errors)}) — pipeline will NOT start:\n")
        for i, c in enumerate(errors, 1):
            print(f"  [{i:02d}] [{c.category}] {c.label}")
            for ln in c.message.splitlines():
                print(f"        {ln}")
            print()

    if warnings:
        print(f"\n  ⚠  WARNINGS  ({len(warnings)}) — pipeline can run but review:\n")
        for c in warnings:
            print(f"  [W]  [{c.category}] {c.label}")
            for ln in c.message.splitlines():
                print(f"        {ln}")
            print()

    print(f"  ✔ PASSED : {len(passed)}/{len(all_checks)} checks")
    if errors:
        print(f"  ✗ FAILED : {len(errors)} error(s)  (BLOCKING)")
    if warnings:
        print(f"  ⚠ WARN   : {len(warnings)} warning(s) (non-blocking)")
    print(f"\n{sep}\n")

    if errors:
        logging.error(
            f"[PREFLIGHT] {len(errors)} blocking error(s) found. "
            f"Fix the issues above and re-run. Pipeline aborted.")
        return 1

    if warnings:
        logging.warning(
            f"[PREFLIGHT] {len(warnings)} warning(s). "
            f"Review them above; continuing...")

    logging.info(f"[PREFLIGHT] All {len(passed)} checks passed. Proceeding.")
    return 0


# ── Tools-only preflight (used by `check` subcommand) ────────────

def preflight_tools_only_check(
        stages: Set[str],
        params: dict,
        out_root: Path,
        use_smk: bool = False) -> int:
    """
    Run only environment / tool checks — no samples, no lineage required.
    Used exclusively by the `check` subcommand.

    Checks performed:
      A. Python version
      B. Conda / mamba availability
      C. Tool binaries (fastp / spades / busco) — version probe
      D. Snakemake (if --smk)
      E. System resources (CPUs, RAM, disk) — informational

    Checks intentionally OMITTED (require real samples / full params):
      F. Parameter sanity (busco_lineage, threads, etc.)
      G. Sample integrity
      H. Output dir writable
      I. BUSCO fasta_dir
      J. Snakemake conda env bindings
    """
    all_checks: List[_Check] = []
    conda_exe = params.get("conda_exe", "conda")

    # A. Python
    all_checks.append(_chk_python())

    # B. Conda (only if any conda_env_* is set)
    needs_conda = any(params.get(k) for k in
                      ("conda_env_fastp", "conda_env_spades", "conda_env_busco"))
    chk_b = _chk_conda(conda_exe, needs_conda)
    if chk_b:
        all_checks.append(chk_b)

    # C. Tool binaries — only for requested stages
    if "fastp" in stages:
        all_checks.append(_probe_tool(
            "fastp", params.get("fastp_exe", "fastp"),
            params.get("conda_env_fastp"), params.get("fastp_path"), conda_exe))
    if "spades" in stages:
        all_checks.append(_probe_tool(
            "spades", params.get("spades_exe", "spades.py"),
            params.get("conda_env_spades"), params.get("spades_path"), conda_exe))
    if "busco" in stages:
        all_checks.append(_probe_tool(
            "busco", params.get("busco_exe", "busco"),
            params.get("conda_env_busco"), params.get("busco_path"), conda_exe))

    # D. Snakemake
    chk_d = _chk_snakemake(use_smk)
    if chk_d:
        all_checks.append(chk_d)

    # E. System resources (informational; disk only if out_root already exists)
    cpus, mem_gb = get_system_resources()
    all_checks.append(_Check("E", "System CPUs/RAM", True,
                             f"CPUs={cpus}  RAM={mem_gb} GB"))
    if out_root.exists():
        free_gb = get_disk_free_gb(out_root)
        if free_gb >= 0:
            disk_ok = free_gb >= 5.0
            all_checks.append(_Check(
                "E", "Disk free space",
                disk_ok,
                f"{free_gb:.1f} GB free on {out_root.anchor}"
                + ("" if disk_ok else " — WARNING: < 5 GB"),
                fatal=False))

    # ── Print report ─────────────────────────────────────────────
    sep = "=" * 66
    print(f"\n{sep}")
    print(f"  ENVIRONMENT CHECK  (v6.25 — {len(all_checks)} checks)")
    print(f"  Stages: {', '.join(sorted(stages)) or 'none'}")
    print(sep)

    errors   = [c for c in all_checks if not c.ok and c.fatal]
    warnings = [c for c in all_checks if not c.ok and not c.fatal]
    passed   = [c for c in all_checks if c.ok]

    if errors:
        print(f"\n  ✗ ERRORS  ({len(errors)}) — fix before running the pipeline:\n")
        for i, c in enumerate(errors, 1):
            print(f"  [{i:02d}] [{c.category}] {c.label}")
            for ln in c.message.splitlines():
                print(f"        {ln}")
            print()

    if warnings:
        print(f"\n  ⚠  WARNINGS  ({len(warnings)}) — non-blocking:\n")
        for c in warnings:
            print(f"  [W]  [{c.category}] {c.label}")
            for ln in c.message.splitlines():
                print(f"        {ln}")
            print()

    for c in passed:
        short_msg = c.message.splitlines()[0] if c.message else ""
        print(f"  ✔  [{c.category}] {c.label}: {short_msg}")

    print(f"\n{sep}")
    print(f"  PASSED  : {len(passed)}/{len(all_checks)}")
    if errors:
        print(f"  ERRORS  : {len(errors)}  ← fix these before running")
    if warnings:
        print(f"  WARNINGS: {len(warnings)}")
    print(f"{sep}\n")

    return 1 if errors else 0


# ── Backward-compat shims (called from old code paths) ───────────

def preflight_env_check(cmd: str, params: dict,
                        strict: bool = True) -> int:
    """Deprecated shim — kept for API compatibility only.
    BUG-FIX 18: previously called preflight_full_check() with fake_samples=[],
    which instantly hit _chk_samples() 'No samples found' -> rc=1 for any
    external caller.  Now delegates to preflight_tools_only_check() exactly
    as the `check` subcommand does, so it works without samples or
    busco_lineage."""
    out_root = Path(params.get("_out_root", ".")).resolve()
    use_smk  = params.get("use_snakemake", False)
    # Determine which stages to probe from the legacy `cmd` argument.
    _cmd_stage_map: Dict[str, Set[str]] = {
        "all":          {"fastp", "spades", "busco"},
        "fastp":        {"fastp"},
        "spades":       {"spades"},
        "busco":        {"busco"},
        "spades2busco": {"spades", "busco"},
    }
    stages: Set[str] = _cmd_stage_map.get(cmd, {"fastp", "spades", "busco"})
    return preflight_tools_only_check(stages, params, out_root, use_smk)


def validate_all_params(cmd: str, params: dict,
                        samples: List["Sample"]) -> List[str]:
    """Deprecated shim — preflight_full_check now does all validation.
    Returns empty list (validation already ran); call site handles rc."""
    return []


# ==============================================================
# 9. SPAdes utilities
# ==============================================================

VALID_ASM_MODES: Set[str] = {
    "isolate", "sc", "plasmid", "meta", "metaplasmid", "metaviral",
    "rna", "rnaviral", "bio", "corona", "sewage",
}


def normalize_assembly_mode(val: str) -> str:
    v = (val or "").strip().strip("'\"")
    return re.sub(r"^-+", "", v)


def spades_expected_outputs(data_type: str) -> List[str]:
    return (["transcripts.fasta", "contigs.fasta", "scaffolds.fasta"]
            if data_type == "rna"
            else ["scaffolds.fasta", "contigs.fasta"])


def spades_success_files(spades_dir: Path, data_type: str) -> List[Path]:
    return [
        p for nm in spades_expected_outputs(data_type)
        for p in [spades_dir / nm]
        if p.exists() and p.stat().st_size > 0
    ]


def spades_has_partial_state(spades_dir: Path) -> bool:
    """Return True ONLY when SPAdes has written a real checkpoint.

    Official SPAdes docs: --continue "starts from the last available
    check-point."  Checkpoints are stored in pipeline_state/ (written
    after each graph-building stage) and K<N>/ subdirectories.

    spades.log / params.txt / dataset.info exist in ALL runs including
    completely failed ones.  Using them as indicators causes --continue
    to be tried on dead runs → rc=255 every time → infinite retry loop.
    NEVER include these files as partial-state indicators.
    """
    if not spades_dir.is_dir():
        return False
    # pipeline_state/ is the authoritative checkpoint written by SPAdes
    # after completing each assembly stage.  Must be non-empty.
    ps = spades_dir / "pipeline_state"
    if ps.is_dir():
        try:
            if any(ps.iterdir()):
                return True
        except PermissionError:
            pass
    # K<N>/ subdirectories mean per-kmer assembly has progressed
    try:
        for p in spades_dir.iterdir():
            if p.is_dir() and re.match(r"^K\d+$", p.name):
                try:
                    if any(p.iterdir()):
                        return True
                except PermissionError:
                    pass
    except PermissionError:
        pass
    return False


def parse_keep_files(raw: Optional[List[str]]) -> List[str]:
    if not raw:
        return []
    items: List[str] = []
    seen: Set[str] = set()
    for v in raw:
        for x in str(v).split(","):
            x = x.strip()
            if x and x not in seen:
                seen.add(x)
                items.append(x)
    return items


def default_keep_patterns() -> List[str]:
    return [
        "contigs.fasta", "scaffolds.fasta", "transcripts.fasta",
        "assembly_graph.fastg", "assembly_graph_with_scaffolds.gfa",
        "assembly_graph_after_simplification.gfa",
        "contigs.paths", "scaffolds.paths", "transcripts.paths",
        "params.txt", "dataset.info", "input_dataset.yaml",
        "run_spades.sh", "run_spades.yaml", "spades.log", "warnings.log",
        "cmd_spades.txt", "stdout_spades.log", "stderr_spades.log",
        "returncode_spades.txt",
        ".spades.done", ".spades.fail",
        "*.log", "*.txt", "*.yaml", "*.yml",
    ]


def spades_clean_dir(spades_dir: Path, keep_patterns: List[str]) -> None:
    for dname in ["tmp", "misc", "logs", "pipeline_state",
                  "corrected", ".bin_reads"]:
        d = spades_dir / dname
        if d.is_dir():
            shutil.rmtree(d, ignore_errors=True)
    for d in list(spades_dir.iterdir()):
        if d.is_dir() and re.match(r"^K\d+$", d.name):
            shutil.rmtree(d, ignore_errors=True)
    keep_set: Set[Path] = set()
    for pat in keep_patterns:
        for p in spades_dir.glob(pat):
            if p.is_file():
                keep_set.add(p.resolve())
    for p in list(spades_dir.iterdir()):
        if p.is_file() and p.resolve() not in keep_set:
            try:
                p.unlink()
            except Exception:
                pass


def spades_wipe_partial(spades_dir: Path) -> None:
    """Wipe SPAdes output dir for a clean restart.

    v6.7: wipe EVERYTHING including spades.log, params.txt, dataset.info.
    Previous versions preserved spades.log which caused spades_has_partial_state()
    to still return True on the next Snakemake retry → infinite --continue loop.

    Our own pipeline stamps (.spades.fail) are kept so we can track failures.
    Previous logs are renamed to PREV_* so diagnostics are not lost.
    """
    if not spades_dir.is_dir():
        return
    # Rename previous diagnostic logs (keep for diagnosis)
    # Use replace() not rename() — replace() overwrites existing PREV_* files
    for log_name in ("stdout_spades.log", "stderr_spades.log",
                     "returncode_spades.txt", "cmd_spades.txt"):
        p = spades_dir / log_name
        if p.exists():
            try:
                p.replace(spades_dir / f"PREV_{log_name}")
            except Exception:
                pass
    # Remove done stamp (if somehow present)
    (spades_dir / ".spades.done").unlink(missing_ok=True)
    # Wipe everything: K<N>/, pipeline_state/, tmp/, misc/, corrected/,
    # .bin_reads/, spades.log, params.txt, dataset.info, assembly_graph*, etc.
    # Keep only: .spades.fail and PREV_* files
    _KEEP = {".spades.fail"}
    for child in list(spades_dir.iterdir()):
        if child.name in _KEEP or child.name.startswith("PREV_"):
            continue
        try:
            if child.is_dir():
                shutil.rmtree(child, ignore_errors=True)
            else:
                child.unlink(missing_ok=True)
        except Exception:
            pass
    logging.warning(
        f"[SPADES_WIPE] {spades_dir.name}: full wipe (spades.log removed "
        f"to prevent false --continue detection; old logs renamed PREV_*)")


# ==============================================================
# 10. Per-sample stage runners
# ==============================================================

def run_fastp_one(sample: Sample, out_root: Path,
                  params: dict, dryrun: bool = False) -> Tuple[bool, str]:
    stage_dir = out_root / params["fastp_dir"] / sample.sid
    stage_dir.mkdir(parents=True, exist_ok=True)
    done = stage_dir / ".fastp.done"
    fail = stage_dir / ".fastp.fail"
    clean_r1 = stage_dir / f"{sample.sid}.R1.clean.fastq.gz"
    clean_r2 = stage_dir / f"{sample.sid}.R2.clean.fastq.gz"

    if fastp_output_complete(stage_dir, sample.sid, done):
        logging.info(f"[FASTP_SKIP] {sample.sid}  (verified complete)")
        return True, "skip"
    if done_ok(done):
        _invalidate_fastp_output(
            stage_dir, sample.sid, done, fail,
            "done_stamp_present_but_json_invalid")
    if (clean_r1.exists() or clean_r2.exists()) and not done_ok(done):
        logging.warning(
            f"[FASTP_STALE] {sample.sid}: output without done stamp -- removing")
        # FIX v6.2: guard with `not dryrun` — previous code deleted real files
        # even in --dryrun mode, making dry-runs destructive.
        if not dryrun:
            for p in [clean_r1, clean_r2,
                      stage_dir / "fastp.json", stage_dir / "fastp.html"]:
                if p.exists():
                    p.unlink(missing_ok=True)

    if not Path(sample.r1).exists():
        msg = f"R1 file not found: {sample.r1}"
        write_fail(fail, msg)
        return False, msg
    if not Path(sample.r2).exists():
        msg = f"R2 file not found: {sample.r2}"
        write_fail(fail, msg)
        return False, msg

    base_args = [
        "-i", sample.r1, "-I", sample.r2,
        "-o", str(clean_r1), "-O", str(clean_r2),
        "-w", str(int(params["threads_fastp"])),
        "-h", str(stage_dir / "fastp.html"),
        "-j", str(stage_dir / "fastp.json"),
    ]
    if int(params.get("fastp_reads_to_process", 0) or 0) > 0:
        base_args += ["--reads_to_process", str(params["fastp_reads_to_process"])]
    if params.get("fastp_detect_adapter_for_pe", True):
        base_args.append("--detect_adapter_for_pe")
    if params.get("data_type") == "rna" and params.get("fastp_trim_poly_x", True):
        base_args.append("--trim_poly_x")
    base_args += extra_args_from_list(params.get("fastp_extra", []))

    try:
        cmd = build_tool_cmd(
            tool_label="FASTP",
            exe_name=params.get("fastp_exe", "fastp"),
            base_args=base_args,
            conda_env=params.get("conda_env_fastp"),
            tool_path=params.get("fastp_path"),
            conda_exe=params.get("conda_exe", "conda"))
    except FileNotFoundError as e:
        write_fail(fail, str(e))
        return False, str(e)

    (stage_dir / "cmd_fastp.txt").write_text(" ".join(cmd) + "\n", encoding="utf-8")
    logging.info(
        f"[FASTP] {sample.sid}  R1={sample.r1}  R2={sample.r2}  "
        f"threads={params['threads_fastp']}")
    rc, out, err = _run_cmd(cmd, cwd=stage_dir, dryrun=dryrun)
    (stage_dir / "stdout_fastp.log").write_text(out, encoding="utf-8", errors="ignore")
    (stage_dir / "stderr_fastp.log").write_text(err, encoding="utf-8", errors="ignore")
    (stage_dir / "returncode_fastp.txt").write_text(str(rc) + "\n", encoding="utf-8")

    if rc != 0:
        write_fail(fail, f"fastp_returncode={rc}")
        done.unlink(missing_ok=True)
        logging.error(f"[FASTP_FAIL] {sample.sid}  rc={rc}")
        return False, f"fastp_returncode={rc}"
    if not clean_r1.exists() or not clean_r2.exists():
        write_fail(fail, "fastp_no_output")
        done.unlink(missing_ok=True)
        return False, "fastp_no_output"

    write_done(done)
    fail.unlink(missing_ok=True)
    return True, "ok"


def run_spades_one(sample: Sample, out_root: Path,
                   params: dict, dryrun: bool = False) -> Tuple[bool, str]:
    spades_root = spades_sample_dir(out_root, params, sample.sid)
    spades_root.mkdir(parents=True, exist_ok=True)
    done = spades_root / ".spades.done"
    fail = spades_root / ".spades.fail"
    data_type = params.get("data_type", "wgs")

    ok_files = spades_success_files(spades_root, data_type)
    if ok_files:
        if not done_ok(done):
            write_done(done)
        logging.info(f"[SPADES_SKIP] {sample.sid}  (outputs exist)")
        return True, "skip"

    if params.get("use_fastp", True):
        fd = out_root / params["fastp_dir"] / sample.sid
        r1 = str((fd / f"{sample.sid}.R1.clean.fastq.gz").resolve())
        r2 = str((fd / f"{sample.sid}.R2.clean.fastq.gz").resolve())
        if not Path(r1).exists():
            msg = f"fastp output missing: {r1}"
            write_fail(fail, msg)
            return False, msg
        if not Path(r2).exists():
            msg = f"fastp output missing: {r2}"
            write_fail(fail, msg)
            return False, msg
    else:
        r1, r2 = sample.r1, sample.r2
        if not Path(r1).exists():
            msg = f"R1 not found: {r1}"
            write_fail(fail, msg)
            return False, msg
        if not Path(r2).exists():
            msg = f"R2 not found: {r2}"
            write_fail(fail, msg)
            return False, msg

    asm_mode = normalize_assembly_mode(params.get("assembly_mode", ""))
    if asm_mode and asm_mode not in VALID_ASM_MODES:
        msg = f"invalid_assembly_mode={asm_mode}"
        write_fail(fail, msg)
        done.unlink(missing_ok=True)
        return False, msg
    asm_flag = f"--{asm_mode}" if asm_mode else ""

    has_partial = spades_has_partial_state(spades_root)
    prev_failed = fail.exists()
    if has_partial and prev_failed:
        logging.warning(
            f"[SPADES_STALE] {sample.sid}: partial from FAILED run -- wiping")
        spades_wipe_partial(spades_root)
        has_partial = False

    def _fresh_args() -> List[str]:
        a = [
            "-1", r1, "-2", r2,
            "-o", str(spades_root),
            "-t", str(int(params["threads_spades"])),
            "-m", str(int(params["spades_mem_gb"])),
            "--only-assembler",
        ]
        if asm_flag:
            a.append(asm_flag)
        elif data_type == "rna":
            a.append("--rna")
        a += extra_args_from_list(params.get("spades_extra", []))
        return a

    def _continue_args() -> List[str]:
        # SPAdes: --continue ONLY accepts -o. No -t, -m, or any other flag.
        # "you cannot specify any option except -o with --continue option!"
        return ["--continue", "-o", str(spades_root)]

    if has_partial:
        base_args = _continue_args()
        logging.info(f"[SPADES_CONTINUE] {sample.sid}  dir={spades_root.name}")
        is_continue = True
    else:
        base_args = _fresh_args()
        logging.info(
            f"[SPADES] {sample.sid}  R1={r1}  R2={r2}  "
            f"threads={params['threads_spades']}  mem={params['spades_mem_gb']}GB"
            f"  dir={spades_root.name}")
        is_continue = False

    try:
        cmd = build_tool_cmd(
            tool_label="SPADES",
            exe_name=params.get("spades_exe", "spades.py"),
            base_args=base_args,
            conda_env=params.get("conda_env_spades"),
            tool_path=params.get("spades_path"),
            conda_exe=params.get("conda_exe", "conda"))
    except FileNotFoundError as e:
        write_fail(fail, str(e))
        return False, str(e)

    def _write_logs(out_txt: str, err_txt: str, rc_val: int) -> None:
        (spades_root / "cmd_spades.txt").write_text(
            " ".join(cmd) + "\n", encoding="utf-8")
        (spades_root / "stdout_spades.log").write_text(
            out_txt, encoding="utf-8", errors="ignore")
        (spades_root / "stderr_spades.log").write_text(
            err_txt, encoding="utf-8", errors="ignore")
        (spades_root / "returncode_spades.txt").write_text(
            str(rc_val) + "\n", encoding="utf-8")

    rc, out, err = _run_cmd(cmd, cwd=spades_root, dryrun=dryrun)
    _write_logs(out, err, rc)

    if rc != 0 and is_continue and not dryrun:
        logging.warning(
            f"[SPADES_CONTINUE_FAIL] {sample.sid}: rc={rc}; retrying fresh")
        spades_wipe_partial(spades_root)
        fresh_args = _fresh_args()
        try:
            cmd = build_tool_cmd(
                tool_label="SPADES",
                exe_name=params.get("spades_exe", "spades.py"),
                base_args=fresh_args,
                conda_env=params.get("conda_env_spades"),
                tool_path=params.get("spades_path"),
                conda_exe=params.get("conda_exe", "conda"))
        except FileNotFoundError as e:
            write_fail(fail, str(e))
            return False, str(e)
        logging.info(f"[SPADES_RETRY] {sample.sid}")
        rc, out, err = _run_cmd(cmd, cwd=spades_root, dryrun=dryrun)
        _write_logs(out, err, rc)

    if rc != 0:
        write_fail(fail, f"spades_returncode={rc}")
        done.unlink(missing_ok=True)
        logging.error(f"[SPADES_FAIL] {sample.sid}  rc={rc}\n{tail_text(err)}")
        return False, f"spades_returncode={rc}"

    if not spades_success_files(spades_root, data_type):
        write_fail(fail, "spades_no_success_output")
        done.unlink(missing_ok=True)
        return False, "spades_no_success_output"

    write_done(done)
    fail.unlink(missing_ok=True)

    if params.get("spades_clean", True) and not dryrun:
        kp = params.get("keep_files") or default_keep_patterns()
        spades_clean_dir(spades_root, kp)
    return True, "ok"


def busco_has_summary(run_dir: Path) -> bool:
    """Return True only when a real short_summary*.txt file exists at the top
    level of *run_dir*.  Using is_file() prevents false positives from
    directories named like short_summary*, and is_dir() avoids calling glob()
    on non-directory paths."""
    return run_dir.is_dir() and any(
        p for p in run_dir.glob("short_summary*.txt") if p.is_file()
    )


def choose_busco_fasta(sample: Sample, out_root: Path,
                       params: dict) -> Path:
    bfd = params.get("busco_fasta_dir")
    if bfd:
        fasta_dir = Path(bfd)
        for ext in (".fasta", ".fa", ".fna"):
            p = fasta_dir / f"{sample.sid}{ext}"
            if p.exists() and p.stat().st_size > 0:
                return p.resolve()
        cands = (list(fasta_dir.glob(f"{sample.sid}*.fasta")) +
                 list(fasta_dir.glob(f"{sample.sid}*.fa")) +
                 list(fasta_dir.glob(f"{sample.sid}*.fna")))
        if cands:
            return cands[0].resolve()
        # --busco_fasta_dir was explicitly given but no matching file found.
        # Return the expected path (which does not exist) so the downstream
        # existence check in run_busco_one() raises a clear error instead of
        # silently falling through to an unrelated SPAdes assembly.
        expected = (fasta_dir / f"{sample.sid}.fasta").resolve()
        logging.warning(
            f"[BUSCO_FASTA] {sample.sid}: --busco_fasta_dir set but no matching "
            f"fasta found in {fasta_dir}. Expected: {expected}")
        return expected
    spades_dir = spades_sample_dir(out_root, params, sample.sid)
    data_type = params.get("data_type", "wgs")
    which = (params.get("which_fasta_to_busco") or "").strip()
    if not which or which.upper() == "AUTO":
        which = "transcripts.fasta" if data_type == "rna" else "scaffolds.fasta"
    wp = Path(which)
    if wp.is_absolute() and wp.is_file():
        return wp.resolve()
    cand = spades_dir / which
    if cand.exists() and cand.stat().st_size > 0:
        return cand.resolve()
    for nm in spades_expected_outputs(data_type):
        q = spades_dir / nm
        if q.exists() and q.stat().st_size > 0:
            return q.resolve()
    return (spades_dir / "scaffolds.fasta").resolve()


def run_busco_one(sample: Sample, out_root: Path,
                  params: dict, dryrun: bool = False) -> Tuple[bool, str]:
    run_dir = busco_run_dir(out_root, params, sample.sid)
    busco_root = run_dir.parent
    run_dir.mkdir(parents=True, exist_ok=True)
    done = run_dir / ".busco.done"
    fail = run_dir / ".busco.fail"

    if busco_has_summary(run_dir) and done_ok(done):
        logging.info(f"[BUSCO_SKIP] {sample.sid}")
        return True, "skip_existing"

    asm = choose_busco_fasta(sample, out_root, params)
    if not asm.exists() or asm.stat().st_size == 0:
        msg = f"missing_assembly_fasta: {asm}"
        write_fail(fail, msg)
        done.unlink(missing_ok=True)
        return False, msg

    data_type = params.get("data_type", "wgs")
    mode = params.get("busco_mode") or \
           ("transcriptome" if data_type == "rna" else "genome")

    busco_ver = detect_busco_version(
        busco_exe=params.get("busco_exe", "busco"),
        conda_env=params.get("conda_env_busco"),
        tool_path=params.get("busco_path"),
        conda_exe=params.get("conda_exe", "conda"),
    )
    major = busco_ver[0]

    is_restart = busco_has_partial_state(run_dir)
    if is_restart:
        logging.info(
            f"[BUSCO_PARTIAL] {sample.sid}: partial run -- attempting --restart")

    base_args = [
        "-i", str(asm),
        "-l", str(params["busco_lineage"]),
        "-m", str(mode),
        "-o", sample.sid,
        "--out_path", str(busco_root),
        "-c", str(int(params["threads_busco"])),
    ]

    if is_restart and not params.get("busco_force", False):
        base_args += _busco_build_restart_args(busco_ver, run_dir)
    if params.get("busco_offline", False):
        base_args.append("--offline")
    if params.get("busco_force", False):
        base_args.append("-f")
    base_args += extra_args_from_list(params.get("busco_extra", []))

    try:
        cmd = build_tool_cmd(
            tool_label="BUSCO",
            exe_name=params.get("busco_exe", "busco"),
            base_args=base_args,
            conda_env=params.get("conda_env_busco"),
            tool_path=params.get("busco_path"),
            conda_exe=params.get("conda_exe", "conda"))
    except FileNotFoundError as e:
        write_fail(fail, str(e))
        return False, str(e)

    action = "BUSCO_RESTART" if is_restart else "BUSCO"
    logging.info(
        f"[{action}] {sample.sid}  fasta={asm.name}  "
        f"v{major}.{busco_ver[1]}  mode={mode}  lineage={params['busco_lineage']}"
        f"  offline={params.get('busco_offline', False)}")
    # BUSCO refuses to start when -o <run_dir> already exists (without -f).
    # On a fresh run run_dir only holds our own stamps/logs, so remove it
    # and let BUSCO create it; logs are written back afterwards.
    if not is_restart and not dryrun:
        shutil.rmtree(run_dir, ignore_errors=True)
    rc, out, err = _run_cmd(cmd, cwd=busco_root, dryrun=dryrun)
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "cmd_busco.txt").write_text(" ".join(cmd) + "\n", encoding="utf-8")
    (run_dir / "stdout_busco.log").write_text(out, encoding="utf-8", errors="ignore")
    (run_dir / "stderr_busco.log").write_text(err, encoding="utf-8", errors="ignore")
    (run_dir / "returncode_busco.txt").write_text(str(rc) + "\n", encoding="utf-8")

    # v6.7: if -r (restart) failed, wipe entirely and retry fresh.
    # Per official BUSCO docs: -f means "force overwrite of existing files",
    # NOT a resume option. Never auto-add -f after -r fails.
    if rc != 0 and is_restart and not dryrun:
        logging.warning(
            f"[BUSCO_RESTART_FAIL] {sample.sid}: -r (restart) rc={rc}. "
            f"Wiping run dir and retrying fresh (NOT -f).")
        try:
            shutil.rmtree(run_dir)
            run_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            logging.warning(f"[BUSCO_WIPE] {sample.sid}: wipe error: {e}")
        fresh_base = [a for a in base_args if a not in ("--restart", "-r")]
        try:
            fresh_cmd_list = build_tool_cmd(
                tool_label="BUSCO",
                exe_name=params.get("busco_exe", "busco"),
                base_args=fresh_base,
                conda_env=params.get("conda_env_busco"),
                tool_path=params.get("busco_path"),
                conda_exe=params.get("conda_exe", "conda"))
        except FileNotFoundError as e:
            write_fail(fail, str(e))
            return False, str(e)
        (run_dir / "cmd_busco.txt").write_text(
            " ".join(fresh_cmd_list) + "\n", encoding="utf-8")
        logging.info(f"[BUSCO_FRESH] {sample.sid}: fresh run after -r failure")
        rc, out, err = _run_cmd(fresh_cmd_list, cwd=busco_root, dryrun=dryrun)
        (run_dir / "stdout_busco.log").write_text(
            out, encoding="utf-8", errors="ignore")
        (run_dir / "stderr_busco.log").write_text(
            err, encoding="utf-8", errors="ignore")
        (run_dir / "returncode_busco.txt").write_text(
            str(rc) + "\n", encoding="utf-8")

    if rc != 0:
        write_fail(fail, f"busco_returncode={rc}")
        done.unlink(missing_ok=True)
        logging.error(f"[BUSCO_FAIL] {sample.sid}  rc={rc}")
        return False, f"busco_returncode={rc}"
    if not busco_has_summary(run_dir):
        write_fail(fail, "busco_no_summary")
        done.unlink(missing_ok=True)
        return False, "busco_no_summary"

    write_done(done)
    fail.unlink(missing_ok=True)
    return True, "ok"


# ==============================================================
# 11. Summary TSV
# ==============================================================

def write_summary_tsv(samples: List[Sample], out_root: Path,
                      params: dict) -> Path:
    summary_path = out_root / params["summary_tsv_name"]
    header = ["sample", "r1", "r2",
              "fastp_done", "spades_done", "busco_done",
              "spades_subdir", "busco_subdir"]
    rows: List[List[str]] = []
    for s in samples:
        fd = "1" if done_ok(
            out_root / params["fastp_dir"] / s.sid / ".fastp.done") else "0"
        sd_dir = spades_sample_dir(out_root, params, s.sid)
        sd = "1" if done_ok(sd_dir / ".spades.done") else "0"
        bd_dir = busco_run_dir(out_root, params, s.sid)
        bd = "1" if done_ok(bd_dir / ".busco.done") else "0"
        rows.append([s.sid, s.r1, s.r2, fd, sd, bd,
                     str(sd_dir.name), str(bd_dir.name)])
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    logging.info(f"[SUMMARY] {summary_path}  rows={len(rows)}")
    return summary_path


# ==============================================================
# 12. Parallel runner (original ThreadPoolExecutor backend)
# ==============================================================

def _run_stage_parallel(
    *, stage_name: str, samples: List[Sample],
    fn, out_root: Path, params: dict,
    max_workers: int, keep_going: bool, dryrun: bool,
) -> Tuple[List[str], List[str]]:
    ok_list: List[str] = []
    fail_list: List[str] = []
    logging.info(
        f"[STAGE] -- {stage_name} --  "
        f"samples={len(samples)}  workers={max_workers}")
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = {ex.submit(fn, s, out_root, params, dryrun): s.sid
                for s in samples}
        for fut in as_completed(futs):
            sid = futs[fut]
            try:
                ok, reason = fut.result()
            except Exception as e:
                ok, reason = False, f"exception: {e}"
            if ok:
                ok_list.append(sid)
            else:
                fail_list.append(sid)
                logging.error(f"[{stage_name}_FAIL] {sid}: {reason}")
                if not keep_going:
                    # cancel() only prevents *queued* tasks from starting;
                    # tasks already running will finish naturally.  We break
                    # out of the as_completed loop so no further results are
                    # processed, but the executor context manager will still
                    # wait for in-flight tasks to complete before returning.
                    for f in futs:
                        f.cancel()
                    break
    logging.info(
        f"[STAGE] -- {stage_name} done --  "
        f"ok={len(ok_list)}  fail={len(fail_list)}")
    return ok_list, fail_list


def run_direct_pipeline(samples: List[Sample], out_root: Path,
                        params: dict, cmd: str, dryrun: bool) -> int:
    keep_going = bool(params.get("keep_going", True))
    max_parallel = max(1, int(params.get("max_parallel", 8)))
    any_fail = False
    summary_tsv = out_root / params.get("summary_tsv_name", "summary_report.tsv")
    raw_skip = params.get("skip_samples_raw")

    def _skip_set(stage: str) -> Set[str]:
        return load_skip_set(
            raw_vals=raw_skip, base_dir=out_root,
            summary_tsv=summary_tsv if summary_tsv.exists() else None,
            stage=stage)

    def _filter_manual(stage: str, slist: List[Sample]) -> List[Sample]:
        ss = _skip_set(stage)
        if not ss:
            return slist
        kept = [s for s in slist if s.sid not in ss]
        skipped = [s.sid for s in slist if s.sid in ss]
        if skipped:
            logging.info(
                f"[SKIP] {stage}: {len(skipped)} skipped (--skip_samples): "
                f"{', '.join(skipped[:20])}"
                + ("..." if len(skipped) > 20 else ""))
        return kept

    if cmd in ("all", "fastp") and params.get("use_fastp", True):
        active = _filter_manual("fastp", samples)
        if active:
            _, active = prescan_fastp(active, out_root, params)
        if active:
            _, fails = _run_stage_parallel(
                stage_name="FASTP", samples=active, fn=run_fastp_one,
                out_root=out_root, params=params,
                max_workers=max_parallel, keep_going=keep_going, dryrun=dryrun)
            if fails:
                any_fail = True
                if not keep_going:
                    return 1
        else:
            logging.info("[STAGE] FASTP -- all samples done / skipped")

    if cmd in ("all", "spades", "spades2busco"):
        active = _filter_manual("spades", samples)
        if active:
            _, active = prescan_spades(active, out_root, params)
        if active:
            _, fails = _run_stage_parallel(
                stage_name="SPADES", samples=active, fn=run_spades_one,
                out_root=out_root, params=params,
                max_workers=max_parallel, keep_going=keep_going, dryrun=dryrun)
            if fails:
                any_fail = True
                if not keep_going:
                    return 1
        else:
            logging.info("[STAGE] SPADES -- all samples done / skipped")

    if cmd in ("all", "busco", "spades2busco"):
        active = _filter_manual("busco", samples)
        if active:
            _, active = prescan_busco(active, out_root, params)
        if active:
            _, fails = _run_stage_parallel(
                stage_name="BUSCO", samples=active, fn=run_busco_one,
                out_root=out_root, params=params,
                max_workers=max_parallel, keep_going=keep_going, dryrun=dryrun)
            if fails:
                any_fail = True
        else:
            logging.info("[STAGE] BUSCO -- all samples done / skipped")

    if not dryrun:
        write_summary_tsv(samples, out_root, params)
    return 1 if any_fail else 0


# ==============================================================
# 12b. Snakemake backend  (FIXED v6.1)
#
#   Key fixes:
#   1) BUSCO fasta moved from input: to params: -- avoids
#      MissingInputException because Snakemake no longer tries
#      to find a rule producing that file in the DAG.
#   2) SPAdes shell: fresh_cmd defined once at top of shell block
#      so both if/else branches can reference it.
#   3) BUSCO partial-state detection uses proper shell glob
#      (compgen -G) instead of [ -d "...run_*" ].
#   4) Added Snakemake --touch for already-complete samples when
#      mixed v5/v6 results exist.
# ==============================================================

def _smk_tool_prefix(tool_key: str, params: dict) -> str:
    exe_map = {"fastp": "fastp_exe", "spades": "spades_exe", "busco": "busco_exe"}
    path_map = {"fastp": "fastp_path", "spades": "spades_path", "busco": "busco_path"}
    env_map = {"fastp": "conda_env_fastp", "spades": "conda_env_spades",
               "busco": "conda_env_busco"}
    exe_name = params.get(exe_map[tool_key], tool_key)
    conda_env = params.get(env_map[tool_key])
    tool_path = params.get(path_map[tool_key])
    conda_exe = params.get("conda_exe", "conda")
    exe_for_shell = exe_name
    if tool_path:
        p = Path(tool_path)
        exe_for_shell = str((p / exe_name) if p.is_dir() else p)
    if conda_env:
        return f"{shlex.quote(str(conda_exe))} run -n {shlex.quote(str(conda_env))} --no-capture-output {shlex.quote(str(exe_for_shell))}"
    return shlex.quote(str(exe_for_shell)) if tool_path else shlex.quote(str(exe_name))


def _assert_no_conda_run_double_dash(cmd: str, label: str = "command") -> str:
    if " run -- " in f" {cmd} " or " run  -- " in f" {cmd} ":
        raise ValueError(
            f"[SNAKEMAKE] Invalid {label}: found forbidden `conda run -- ...` form.\n"
            f"  Command: {cmd}\n"
            f"  Fix: remove the literal `--` after `conda run`."
        )
    return cmd


def _auto_smk_spades_parallel(params: dict, cmd: str) -> int:
    """Conservative SPAdes parallelism for embedded Snakemake.

    Rationale:
      - SPAdes is much heavier than fastp/BUSCO in RAM, I/O, temp usage.
      - In user reports, multiple concurrent SPAdes jobs under Snakemake often
        appear to "hang" and then fail in bursts with rc=2.
      - For robustness, cap SPAdes concurrency automatically unless the user
        explicitly lowers/raises max_parallel.
    """
    run_spades = cmd in ("all", "spades", "spades2busco")
    user_mp = max(1, int(params.get("max_parallel", 8) or 8))
    sys_cpus = os.cpu_count() or 1
    thr = max(1, int(params.get("threads_spades", 1) or 1))
    mem_gb = max(1, int(params.get("spades_mem_gb", 1) or 1))

    cpu_cap = max(1, sys_cpus // thr)

    mem_cap = 1
    try:
        with open("/proc/meminfo", encoding="utf-8") as f:
            for ln in f:
                if ln.startswith("MemTotal"):
                    total_mb = int(ln.split()[1]) // 1024
                    safe_mb = max(1024, int(total_mb * 0.80))
                    mem_cap = max(1, safe_mb // max(1024, mem_gb * 1024))
                    break
    except Exception:
        mem_cap = max(1, user_mp)

    if run_spades:
        # Robust default: serialize SPAdes jobs in embedded mode.
        # This is slower, but much more stable for heavy multi-env WGS runs.
        hard_cap = 1
        return max(1, min(user_mp, cpu_cap, mem_cap, hard_cap))
    return max(1, min(user_mp, cpu_cap))


def _smk_extra_fastp(params: dict) -> str:
    parts: List[str] = []
    if int(params.get("fastp_reads_to_process", 0) or 0) > 0:
        parts += ["--reads_to_process", str(params["fastp_reads_to_process"])]
    if params.get("fastp_detect_adapter_for_pe", True):
        parts.append("--detect_adapter_for_pe")
    if params.get("data_type") == "rna" and params.get("fastp_trim_poly_x", True):
        parts.append("--trim_poly_x")
    parts += extra_args_from_list(params.get("fastp_extra", []))
    return " ".join(parts)


def _smk_extra_spades(params: dict) -> str:
    parts: List[str] = ["--only-assembler"]
    asm_mode = normalize_assembly_mode(params.get("assembly_mode", ""))
    if asm_mode and asm_mode in VALID_ASM_MODES:
        parts.append(f"--{asm_mode}")
    elif params.get("data_type") == "rna":
        parts.append("--rna")
    parts += extra_args_from_list(params.get("spades_extra", []))
    return " ".join(parts)


def _smk_extra_busco(params: dict) -> str:
    parts: List[str] = []
    if params.get("busco_offline", False):
        parts.append("--offline")
    if params.get("busco_force", False):
        parts.append("-f")
    parts += extra_args_from_list(params.get("busco_extra", []))
    return " ".join(parts)


def _smk_busco_fasta_name(params: dict) -> str:
    data_type = params.get("data_type", "wgs")
    which = (params.get("which_fasta_to_busco") or "").strip()
    if not which or which.upper() == "AUTO":
        return "transcripts.fasta" if data_type == "rna" else "scaffolds.fasta"
    return which


def generate_snakefile(samples: List[Sample], out_root: Path,
                       params: dict, cmd: str,
                       status_rows: Optional[List[Dict[str, str]]] = None) -> str:
    """Return the text of a self-contained Snakefile.

    v6.10: complete rewrite — use template strings instead of L.append()
    to eliminate all brace-escaping bugs once and for all.
    """
    out          = str(out_root)
    fastp_dir    = params["fastp_dir"]
    spades_dir   = params["spades_dir"]
    busco_dir    = params["busco_dir"]
    suffix       = params.get("suffix_spades_out", "") or ""
    use_fastp    = params.get("use_fastp", True) and cmd in ("all", "fastp")
    run_fastp    = cmd in ("all", "fastp") and use_fastp
    run_spades   = cmd in ("all", "spades", "spades2busco")
    run_busco    = cmd in ("all", "busco",  "spades2busco")
    busco_fasta_dir  = params.get("busco_fasta_dir") or ""
    data_type    = params.get("data_type", "wgs")

    threads_fastp  = int(params["threads_fastp"])
    threads_spades = int(params["threads_spades"])
    threads_busco  = int(params["threads_busco"])
    spades_mem     = int(params["spades_mem_gb"])

    mem_mb_fastp  = int(params.get("smk_mem_mb_fastp",  8_000))
    mem_mb_spades = int(params.get("smk_mem_mb_spades") or (spades_mem * 1024))
    mem_mb_busco  = int(params.get("smk_mem_mb_busco",  16_000))
    disk_mb       = int(params.get("smk_disk_mb",       200_000))

    fastp_prefix  = _smk_tool_prefix("fastp",  params)
    spades_prefix = _smk_tool_prefix("spades", params)
    busco_prefix  = _smk_tool_prefix("busco",  params)

    extra_fastp  = _smk_extra_fastp(params)
    extra_spades = _smk_extra_spades(params)
    extra_busco  = _smk_extra_busco(params)

    busco_lineage    = params.get("busco_lineage", "")
    busco_mode       = params.get("busco_mode") or (
        "transcriptome" if data_type == "rna" else "genome")
    busco_fasta_name = _smk_busco_fasta_name(params)

    # ── SAMPLES dict ─────────────────────────────────────────────────
    samples_dict: Dict[str, Dict[str, str]] = {
        s.sid: {"r1": s.r1, "r2": s.r2} for s in samples
    }

    # ── fastp command ────────────────────────────────────────────────
    # These strings are inserted directly into the Snakefile text via
    # {fastp_cmd} in the f-string template below.
    # Use single braces here so the Snakefile contains {input.r1} etc.
    # which Snakemake then expands at runtime.
    fastp_cmd = (f'{fastp_prefix} '
                 '-i {input.r1} -I {input.r2} '
                 '-o {output.r1_out} -O {output.r2_out} '
                 '-w {threads} '
                 '-h {output.html_f} -j {output.json_f}')
    if extra_fastp:
        fastp_cmd += ' ' + extra_fastp
    fastp_cmd = _assert_no_conda_run_double_dash(fastp_cmd, "fastp_cmd")

    # ── spades commands ───────────────────────────────────────────────
    # FRESH: all args; CONT: ONLY -o (SPAdes v4 rule — no other flags)
    fresh_cmd = (f'{spades_prefix} '
                 '-1 {input.r1} -2 {input.r2} '
                 '-o $SDIR -t {threads} -m {params.mem} '
                 + extra_spades)
    cont_cmd  = f'{spades_prefix} --continue -o $SDIR'
    fresh_cmd = _assert_no_conda_run_double_dash(fresh_cmd, "spades_fresh_cmd")
    cont_cmd  = _assert_no_conda_run_double_dash(cont_cmd,  "spades_continue_cmd")

    # ── busco command ─────────────────────────────────────────────────
    busco_cmd = (f'{busco_prefix} '
                 '-i {params.fasta} '
                 f'-l {busco_lineage} '
                 f'-m {busco_mode} '
                 '-o {wildcards.sample} '
                 '--out_path {params.busco_root} '
                 '-c {threads}')
    if extra_busco:
        busco_cmd += ' ' + extra_busco
    busco_cmd = _assert_no_conda_run_double_dash(busco_cmd, "busco_cmd")

    # ── spades wipe logic (inlined — no bash functions, no brace issues) ──
    # CRITICAL: Do NOT mv/rename stdout_spades.log or stderr_spades.log here!
    # Snakemake opens {log.out} and {log.err} BEFORE the shell block runs and
    # holds file descriptors to these inodes.  Renaming them causes Snakemake's
    # fds to point at the old (now PREV_*) inodes while the subsequent redirect
    # creates new files — Snakemake then sees empty logs and may report failure.
    # The redirect '> {log.out}' will truncate and overwrite the log files anyway.
    spades_wipe_inline = """\
        # Wipe SPAdes internals WITHOUT deleting $SDIR or Snakemake-managed logs
        for _s in pipeline_state tmp misc corrected logs .bin_reads; do
            [ -d "$SDIR/$_s" ] && rm -rf "$SDIR/$_s" 2>/dev/null || true
        done
        for _k in "$SDIR"/K[0-9]*; do
            [ -d "$_k" ] && rm -rf "$_k" 2>/dev/null || true
        done
        rm -f "$SDIR"/spades.log "$SDIR"/params.txt "$SDIR"/dataset.info \\
              "$SDIR"/input_dataset.yaml "$SDIR"/run_spades.sh "$SDIR"/run_spades.yaml \\
              "$SDIR"/assembly_graph* "$SDIR"/scaffolds.* "$SDIR"/contigs.* \\
              "$SDIR"/transcripts.* "$SDIR"/warnings.log "$SDIR"/before_rr.fasta \\
              "$SDIR"/first_pe_contigs.fasta "$SDIR"/returncode_spades.txt \\
              "$SDIR"/cmd_spades.txt "$SDIR"/.spades.done 2>/dev/null || true"""

    # ── busco wipe logic (inlined) ────────────────────────────────────
    # CRITICAL: Same as spades — do NOT mv Snakemake-managed log files.
    busco_wipe_inline = """\
        # Wipe BUSCO internals WITHOUT deleting $RDIR or Snakemake-managed logs
        for _s in "$RDIR"/run_*/ "$RDIR"/tmp_*/; do
            [ -d "$_s" ] && rm -rf "$_s" 2>/dev/null || true
        done
        rm -f "$RDIR"/short_summary*.txt "$RDIR"/*.json "$RDIR"/cmd_busco.txt \\
              "$RDIR"/returncode_busco.txt "$RDIR"/.busco.done 2>/dev/null || true"""
    if data_type == "rna":
        asm_check = ('        if [ ! -s "$SDIR/transcripts.fasta" ] '
                     '&& [ ! -s "$SDIR/contigs.fasta" ]; then')
    else:
        asm_check = ('        if [ ! -s "$SDIR/scaffolds.fasta" ] '
                     '&& [ ! -s "$SDIR/contigs.fasta" ]; then')

    # ── _spades_r1/r2 helper body ─────────────────────────────────────
    if use_fastp and run_spades:
        spades_r1_body = '    return os.path.join(OUT, FASTP_DIR, wc.sample, wc.sample + ".R1.clean.fastq.gz")'
        spades_r2_body = '    return os.path.join(OUT, FASTP_DIR, wc.sample, wc.sample + ".R2.clean.fastq.gz")'
    else:
        spades_r1_body = '    return SAMPLES[wc.sample]["r1"]'
        spades_r2_body = '    return SAMPLES[wc.sample]["r2"]'

    # ── _busco_fasta_path helper body ─────────────────────────────────
    if busco_fasta_dir:
        bfp_lines = [
            f'    bd = {busco_fasta_dir!r}',
            '    for ext in (".fasta", ".fa", ".fna"):',
            '        p = os.path.join(bd, wc.sample + ext)',
            '        if os.path.isfile(p) and os.path.getsize(p) > 0: return p',
            '    cands = sorted(glob.glob(os.path.join(bd, wc.sample + "*")))',
            '    if cands: return cands[0]',
            '    # --busco_fasta_dir set but file not found: return expected path',
            '    # so the shell [ ! -s ... ] guard will surface a clear error.',
            '    return os.path.join(bd, wc.sample + ".fasta")',
        ]
    else:
        fallbacks = (["transcripts.fasta", "contigs.fasta", "scaffolds.fasta"]
                     if data_type == "rna"
                     else ["scaffolds.fasta", "contigs.fasta"])
        bfp_lines = [
            '    spades_d = os.path.join(OUT, SPADES_DIR, wc.sample + SUFFIX)',
            '    primary  = os.path.join(spades_d, BUSCO_FASTA_NAME)',
            '    if os.path.isfile(primary) and os.path.getsize(primary) > 0: return primary',
        ]
        for fb in fallbacks:
            bfp_lines += [
                f'    fb = os.path.join(spades_d, {fb!r})',
                '    if os.path.isfile(fb) and os.path.getsize(fb) > 0: return fb',
            ]
        bfp_lines.append('    return primary')
    busco_fasta_path_body = '\n'.join(bfp_lines)

    # ── rule all target ───────────────────────────────────────────────
    if run_busco:
        rule_all_input = ('        expand(os.path.join(OUT, BUSCO_DIR, "{sample}", ".busco.done"),'
                          ' sample=SAMPLE_IDS),')
    elif run_spades:
        rule_all_input = ('        expand(os.path.join(OUT, SPADES_DIR, "{sample}" + SUFFIX, ".spades.done"),'
                          ' sample=SAMPLE_IDS),')
    elif run_fastp:
        rule_all_input = ('        expand(os.path.join(OUT, FASTP_DIR, "{sample}", ".fastp.done"),'
                          ' sample=SAMPLE_IDS),')
    else:
        rule_all_input = '        [],'

    # ── spades input in rule (depends on whether spades_done needed) ──
    busco_spades_dep = ''
    if run_spades and run_busco and not busco_fasta_dir:
        busco_spades_dep = ('        spades_done = os.path.join(OUT, SPADES_DIR,'
                            ' "{sample}" + SUFFIX, ".spades.done"),')

    # ══════════════════════════════════════════════════════════════════
    # Build Snakefile text section by section
    # NOTE: In this template, {{ and }} are literal { } in the output.
    #       Snakemake placeholders like {wildcards.sample} are written
    #       as plain strings here because we are building the Snakefile
    #       text, NOT running inside a Snakemake shell block.
    # ══════════════════════════════════════════════════════════════════

    parts: List[str] = []

    # ── header + constants ────────────────────────────────────────────
    parts.append(f"""\
import os, glob

SAMPLES = {_json.dumps(samples_dict, indent=2)}
SAMPLE_IDS = sorted(SAMPLES.keys())

OUT              = {out!r}
FASTP_DIR        = {fastp_dir!r}
SPADES_DIR       = {spades_dir!r}
BUSCO_DIR        = {busco_dir!r}
SUFFIX           = {suffix!r}
BUSCO_FASTA_DIR  = {busco_fasta_dir!r}
BUSCO_FASTA_NAME = {busco_fasta_name!r}

def _spades_r1(wc):
{spades_r1_body}

def _spades_r2(wc):
{spades_r2_body}

def _busco_fasta_path(wc):
    \"\"\"Fasta path for BUSCO — in params to avoid DAG resolution.\"\"\"
{busco_fasta_path_body}

rule all:
    input:
{rule_all_input}

localrules: all

wildcard_constraints:
    sample = "[^/]+"
""")

    # ── rule fastp ────────────────────────────────────────────────────
    if run_fastp:
        # Build the fastp shell line separately (contains {input.r1} etc.
        # which must reach the Snakefile as-is for Snakemake to expand)
        fastp_shell_line = f'        {fastp_cmd} > {{log.out}} 2> {{log.err}}'
        # Now build the rule text using plain string — NO f-string —
        # so Python does NOT try to expand {input.r1}, {output.r1_out} etc.
        parts.append(
            "rule fastp:\n"
            "    input:\n"
            '        r1 = lambda wc: SAMPLES[wc.sample]["r1"],\n'
            '        r2 = lambda wc: SAMPLES[wc.sample]["r2"],\n'
            "    output:\n"
            '        r1_out = os.path.join(OUT, FASTP_DIR, "{sample}", "{sample}.R1.clean.fastq.gz"),\n'
            '        r2_out = os.path.join(OUT, FASTP_DIR, "{sample}", "{sample}.R2.clean.fastq.gz"),\n'
            '        json_f = os.path.join(OUT, FASTP_DIR, "{sample}", "fastp.json"),\n'
            '        html_f = os.path.join(OUT, FASTP_DIR, "{sample}", "fastp.html"),\n'
            '        done   = os.path.join(OUT, FASTP_DIR, "{sample}", ".fastp.done"),\n'
            f"    threads: {threads_fastp}\n"
            "    resources:\n"
            f"        mem_mb  = {mem_mb_fastp},\n"
            f"        disk_mb = {disk_mb},\n"
            "    log:\n"
            '        out = os.path.join(OUT, FASTP_DIR, "{sample}", "stdout_fastp.log"),\n'
            '        err = os.path.join(OUT, FASTP_DIR, "{sample}", "stderr_fastp.log"),\n'
            "    params:\n"
            '        stage_dir = os.path.join(OUT, FASTP_DIR, "{sample}"),\n'
            "    shell:\n"
            '        """\n'
            '        set -euo pipefail\n'
            '        echo "================================================================="\n'
            '        echo "[FASTP] sample   : {wildcards.sample}"\n'
            '        echo "[FASTP] R1 input : {input.r1}"\n'
            '        echo "[FASTP] R2 input : {input.r2}"\n'
            '        echo "[FASTP] threads  : {threads}"\n'
            '        echo "[FASTP] mem_mb   : {resources.mem_mb}"\n'
            '        echo "================================================================="\n'
            "        set +e\n"
            + fastp_shell_line + "\n"
            "        FASTP_RC=$?\n"
            "        set -e\n"
            '        echo "$FASTP_RC" > "{params.stage_dir}/returncode_fastp.txt"\n'
            "        if [ $FASTP_RC -ne 0 ]; then exit $FASTP_RC; fi\n"
            '        echo "[FASTP] {wildcards.sample} DONE"\n'
            '        echo "ok=1" > {output.done}\n'
            '        """\n'
            "\n"
        )

    # ── rule spades ───────────────────────────────────────────────────
    if run_spades:
        # Shell lines containing Snakemake placeholders are built outside
        # any f-string so Python does NOT try to expand {input.r1} etc.
        # FIX: fresh_cmd and cont_cmd are injected DIRECTLY as literal
        # shell text — never stored in a Bash variable — to prevent
        # word-splitting on arguments that contain spaces (Bug 3).
        # FIX: Restored checkpoint-detect -> --continue -> fallback logic
        # matching the ThreadPool backend (run_spades_one).  The old
        # 'always-fresh' v6.10 behaviour silently lost the resume path.
        cont_log   = f'            {cont_cmd} > {{log.out}} 2> {{log.err}}'
        fresh_log1 = f'                {fresh_cmd} > {{log.out}} 2> {{log.err}}'
        fresh_log2 = f'            {fresh_cmd} > {{log.out}} 2> {{log.err}}'
        done_line  = '        echo "ok=1" > {output.done}'

        parts.append(
            "rule spades:\n"
            "    input:\n"
            "        r1 = _spades_r1,\n"
            "        r2 = _spades_r2,\n"
            "    output:\n"
            '        done = os.path.join(OUT, SPADES_DIR, "{sample}" + SUFFIX, ".spades.done"),\n'
            f"    threads: {threads_spades}\n"
            "    resources:\n"
            f"        mem_mb  = {mem_mb_spades},\n"
            f"        disk_mb = {disk_mb},\n"
            "    log:\n"
            '        out = os.path.join(OUT, SPADES_DIR, "{sample}" + SUFFIX, "stdout_spades.log"),\n'
            '        err = os.path.join(OUT, SPADES_DIR, "{sample}" + SUFFIX, "stderr_spades.log"),\n'
            "    params:\n"
            '        spades_dir = os.path.join(OUT, SPADES_DIR, "{sample}" + SUFFIX),\n'
            f"        mem        = {spades_mem},\n"
            "    shell:\n"
            '        """\n'
            '        set -euo pipefail\n'
            '        echo "================================================================="\n'
            '        echo "[SPADES] sample  : {wildcards.sample}"\n'
            '        echo "[SPADES] R1      : {input.r1}"\n'
            '        echo "[SPADES] R2      : {input.r2}"\n'
            '        echo "[SPADES] threads : {threads}"\n'
            '        echo "[SPADES] memory  : {params.mem} GB"\n'
            '        echo "[SPADES] mem_mb  : {resources.mem_mb}"\n'
            '        echo "[SPADES] outdir  : {params.spades_dir}"\n'
            '        echo "================================================================="\n'
            '        SDIR="{params.spades_dir}"\n'
            '        mkdir -p "$SDIR"\n'
            "\n"
            "        set +e\n"
            '        if [ -d "$SDIR/pipeline_state" ] && [ -n "$(ls -A "$SDIR/pipeline_state" 2>/dev/null)" ]; then\n'
            '            echo "[SPADES] Checkpoint found -> trying --continue" >&2\n'
            + cont_log + "\n"
            "            SPADES_RC=$?\n"
            "            if [ $SPADES_RC -ne 0 ]; then\n"
            '                echo "[SPADES] --continue failed (rc=$SPADES_RC). Wiping and retrying fresh." >&2\n'
            + spades_wipe_inline + "\n"
            + fresh_log1 + "\n"
            "                SPADES_RC=$?\n"
            "            fi\n"
            "        else\n"
            + spades_wipe_inline + "\n"
            + fresh_log2 + "\n"
            "            SPADES_RC=$?\n"
            "        fi\n"
            "        set -e\n"
            '        echo "$SPADES_RC" > "$SDIR/returncode_spades.txt"\n'
            "        if [ $SPADES_RC -ne 0 ]; then\n"
            '            echo "[SPADES] ERROR: spades exited rc=$SPADES_RC" >&2\n'
            '            if [ -s {log.err} ]; then echo "[SPADES] stderr tail (final failure):" >&2; tail -n 80 {log.err} >&2 || true; fi\n'
            "            exit $SPADES_RC\n"
            "        fi\n"
            + asm_check + "\n"
            '            echo "[SPADES] ERROR: no assembly output" >&2\n'
            "            exit 1\n"
            "        fi\n"
            '        echo "[SPADES] {wildcards.sample} DONE"\n'
            + done_line + "\n"
            '        """\n'
            "\n"
        )

    # ── rule busco ────────────────────────────────────────────────────
    if run_busco:
        busco_r_log    = f'            {busco_cmd} -r > {{log.out}} 2> {{log.err}}'
        busco_fresh1   = f'                {busco_cmd} > {{log.out}} 2> {{log.err}}'
        busco_fresh2   = f'            {busco_cmd} > {{log.out}} 2> {{log.err}}'
        busco_done_line = '        echo "ok=1" > {output.done}'

        parts.append(
            "rule busco:\n"
            # Only emit the "input:" block when there is a real dependency.
            # An empty "input:" block is a Snakemake SyntaxError.
            + ("    input:\n" + busco_spades_dep + "\n" if busco_spades_dep else "")
            + "    output:\n"
            '        done = os.path.join(OUT, BUSCO_DIR, "{sample}", ".busco.done"),\n'
            f"    threads: {threads_busco}\n"
            "    resources:\n"
            f"        mem_mb  = {mem_mb_busco},\n"
            f"        disk_mb = {disk_mb},\n"
            "    log:\n"
            '        out = os.path.join(OUT, BUSCO_DIR, "{sample}", "stdout_busco.log"),\n'
            '        err = os.path.join(OUT, BUSCO_DIR, "{sample}", "stderr_busco.log"),\n'
            "    params:\n"
            "        fasta      = _busco_fasta_path,\n"
            '        run_dir    = os.path.join(OUT, BUSCO_DIR, "{sample}"),\n'
            f"        busco_root = os.path.join(OUT, BUSCO_DIR),\n"
            "    shell:\n"
            '        """\n'
            "        set -euo pipefail\n"
            '        if [ ! -s "{params.fasta}" ]; then\n'
            '            echo "[BUSCO] ERROR: fasta not found or empty: {params.fasta}" >&2\n'
            "            exit 1\n"
            "        fi\n"
            '        echo "================================================================="\n'
            '        echo "[BUSCO] sample  : {wildcards.sample}"\n'
            '        echo "[BUSCO] fasta   : {params.fasta}"\n'
            f'        echo "[BUSCO] lineage : {busco_lineage}"\n'
            f'        echo "[BUSCO] mode    : {busco_mode}"\n'
            '        echo "[BUSCO] threads : {threads}"\n'
            '        echo "[BUSCO] mem_mb  : {resources.mem_mb}"\n'
            '        echo "================================================================="\n'
            '        RDIR="{params.run_dir}"\n'
            '        mkdir -p "$RDIR"\n'
            "\n"
            "        HAS_BUSCO_CHECKPOINT=false\n"
            "        # POSIX-compatible checkpoint detection (compgen is bash-only,\n"
            "        # not available in dash/sh on Debian/Ubuntu). Use a for-glob\n"
            "        # loop with [ -e ] guard instead, which is POSIX-compliant.\n"
            '        for RD in "$RDIR"/run_*/; do\n'
            '            [ -e "$RD" ] || continue\n'
            '            if [ -d "$RD" ] && [ -n "$(find "$RD" -mindepth 1 -maxdepth 2 \\\n'
            '                    -not -name "*.log" -not -name "logs" 2>/dev/null | head -1)" ]; then\n'
            "                HAS_BUSCO_CHECKPOINT=true\n"
            "                break\n"
            "            fi\n"
            "        done\n"
            "\n"
            "        set +e\n"
            '        if [ "$HAS_BUSCO_CHECKPOINT" = "true" ]; then\n'
            '            echo "[BUSCO] Checkpoint found -> -r (restart)"\n'
            + busco_r_log + "\n"
            "            BUSCO_RC=$?\n"
            "            if [ $BUSCO_RC -ne 0 ]; then\n"
            '                echo "[BUSCO] -r failed (rc=$BUSCO_RC). Wiping + fresh run." >&2\n'
            '                touch "$RDIR/.busco.fail" 2>/dev/null || true\n'
            + busco_wipe_inline + "\n"
            + busco_fresh1 + "\n"
            "                BUSCO_RC=$?\n"
            "            fi\n"
            "        else\n"
            + busco_wipe_inline + "\n"
            + busco_fresh2 + "\n"
            "            BUSCO_RC=$?\n"
            "        fi\n"
            "        set -e\n"
            '        echo "$BUSCO_RC" > "$RDIR/returncode_busco.txt"\n'
            "        if [ $BUSCO_RC -ne 0 ]; then\n"
            '            echo "[BUSCO] ERROR: busco exited rc=$BUSCO_RC" >&2\n'
            "            exit $BUSCO_RC\n"
            "        fi\n"
            '        if ! ls "$RDIR"/short_summary*.txt 1>/dev/null 2>&1; then\n'
            '            echo "[BUSCO] ERROR: no summary file produced" >&2\n'
            "            exit 1\n"
            "        fi\n"
            '        echo "[BUSCO] {wildcards.sample} DONE"\n'
            + busco_done_line + "\n"
            '        """\n'
            "\n"
        )

    return "".join(parts)



def _write_v5_compat_stamps(samples: List[Sample], out_root: Path,
                            params: dict) -> int:
    """Write .done stamps for v5 results that are missing them.

    This ensures Snakemake recognises already-complete samples.
    Returns count of stamps written.
    """
    count = 0
    data_type = params.get("data_type", "wgs")

    for s in samples:
        # SPAdes: has output files but no stamp
        sd = spades_sample_dir(out_root, params, s.sid)
        stamp = sd / ".spades.done"
        if not done_ok(stamp) and spades_success_files(sd, data_type):
            write_done(stamp)
            count += 1
            logging.debug(f"[V5_COMPAT] wrote .spades.done for {s.sid}")

        # BUSCO: has summary but no stamp
        rd = busco_run_dir(out_root, params, s.sid)
        stamp = rd / ".busco.done"
        if not done_ok(stamp) and busco_has_summary(rd):
            write_done(stamp)
            count += 1
            logging.debug(f"[V5_COMPAT] wrote .busco.done for {s.sid}")

        # fastp: has output but no stamp
        fd = out_root / params["fastp_dir"] / s.sid
        stamp = fd / ".fastp.done"
        r1 = fd / f"{s.sid}.R1.clean.fastq.gz"
        r2 = fd / f"{s.sid}.R2.clean.fastq.gz"
        json_p = fd / "fastp.json"
        if not done_ok(stamp) and r1.exists() and r2.exists() and json_p.exists():
            try:
                data = _json.loads(json_p.read_text(encoding="utf-8", errors="ignore"))
                if "summary" in data:
                    write_done(stamp)
                    count += 1
                    logging.debug(f"[V5_COMPAT] wrote .fastp.done for {s.sid}")
            except Exception:
                pass

    if count:
        logging.info(f"[V5_COMPAT] Wrote {count} missing .done stamps for pre-existing results")
    return count


def _smk_is_locked(out_root: Path) -> bool:
    """Return True if a Snakemake lock directory exists for *out_root*."""
    lock_dir = out_root / ".snakemake" / "locks"
    if not lock_dir.is_dir():
        return False
    return any(lock_dir.iterdir())


def _smk_force_unlock_files(out_root: Path) -> int:
    """Delete .snakemake/locks/* directly — no Snakefile needed.
    Works even when the Snakefile has parse errors (v6.7).
    """
    lock_dir = out_root / ".snakemake" / "locks"
    if not lock_dir.is_dir():
        logging.info("[SMK_UNLOCK] No lock directory found.")
        return 0
    removed, failed = 0, []
    for item in list(lock_dir.iterdir()):
        try:
            if item.is_file() or item.is_symlink():
                item.unlink()
            else:
                shutil.rmtree(item)
            removed += 1
        except Exception as e:
            failed.append(f"{item.name}: {e}")
    if failed:
        detail = "\n".join(f"  {f}" for f in failed)
        logging.error(
            f"[SMK_UNLOCK] Could not remove {len(failed)} lock file(s):\n"
            f"{detail}\n"
            f"  Manual fix: rm -rf {lock_dir}/*")
        return 1
    logging.info(f"[SMK_UNLOCK] Removed {removed} lock file(s) from {lock_dir}")
    return 0


def _smk_unlock(smk_exe: str, snakefile_path: Path,
                out_root: Path) -> int:
    """Unlock stale Snakemake directory (two-step, v6.7).

    Step 1 (primary): delete lock files directly via Python.
      Works even if the Snakefile has parse errors.
    Step 2 (verify): run snakemake --unlock for bookkeeping.
      If this fails, log a warning only — files already removed.
    """
    rc = _smk_force_unlock_files(out_root)
    if rc != 0:
        return rc
    unlock_cmd = [smk_exe, "--snakefile", str(snakefile_path), "--unlock"]
    logging.info(f"[SMK_UNLOCK] Verifying: {' '.join(unlock_cmd)}")
    try:
        proc = subprocess.run(
            unlock_cmd, cwd=str(out_root),
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, timeout=30)
        for line in (proc.stdout or "").splitlines():
            if line.strip():
                logging.debug(f"[SMK_UNLOCK] {line}")
        if proc.returncode == 0:
            logging.info("[SMK_UNLOCK] Lock cleared (snakemake confirmed).")
        else:
            logging.warning(
                f"[SMK_UNLOCK] snakemake --unlock rc={proc.returncode} "
                f"(lock files already deleted — non-fatal).")
    except Exception as e:
        logging.warning(f"[SMK_UNLOCK] snakemake --unlock skipped: {e}")
    return 0


def _smk_run(smk_cmd: List[str], out_root: Path,
             log_path: Path) -> int:
    """Stream-execute a snakemake command, tee to *log_path*. Return rc."""
    cmd_str = " ".join(smk_cmd)
    logging.info(f"[SMK] Running: {cmd_str}")
    with open(str(log_path), "a", encoding="utf-8") as log_fh:
        log_fh.write(f"CMD: {cmd_str}\n")
        proc = subprocess.Popen(
            smk_cmd, cwd=str(out_root),
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if proc.stdout is None:
            raise RuntimeError(
                "[SMK] subprocess Popen returned stdout=None unexpectedly. "
                "This should not happen with stdout=PIPE.")
        for line in proc.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log_fh.write(line)
        proc.wait()
    return proc.returncode


def run_snakemake_pipeline(samples: List[Sample], out_root: Path,
                           params: dict, cmd: str, dryrun: bool) -> int:
    """Generate a Snakefile and invoke snakemake.

    v6.5 improvements over v6.4:
    - Full status scan before Snakefile generation:
        * Print per-sample status table (DONE/PARTIAL/TODO)
        * Exclude fully-done samples from SAMPLE_IDS (fast)
        * Pass priority map so partial jobs run first
    - resources/priority/wildcard_constraints/localrules in generated rules
    - Auto-unlock and --smk_extra passthrough retained from v6.4
    """

    # 1. Check snakemake availability
    smk_exe = which_or_none("snakemake")
    if not smk_exe:
        logging.error(
            "[SMK] snakemake not found in PATH.\n"
            "  Install: pip install snakemake   or   conda install -c bioconda snakemake\n"
            "  Or run without --smk to use the built-in parallel backend.")
        return 1

    # 1b. Write .done stamps for v5-era results
    _write_v5_compat_stamps(samples, out_root, params)

    # ── v6.5: Full status scan ────────────────────────────────────
    # Classify every sample → DONE / RESUME / TODO
    # Status table already printed by main() — only scan here for
    # filtering and priority logic.
    status_rows = scan_all_status(samples, out_root, params, cmd)

    pending_sids = {r["sid"] for r in status_rows if r["overall"] != _ST_DONE}
    pending_samples = [s for s in samples if s.sid in pending_sids]
    n_skip = len(samples) - len(pending_samples)
    if n_skip:
        logging.info(f"[SMK] {n_skip} sample(s) already complete — excluded from Snakefile")
    if not pending_samples:
        logging.info("[SMK] All samples are complete. Nothing to do.")
        if not dryrun:
            write_summary_tsv(samples, out_root, params)
        return 0
    logging.info(f"[SMK] Pending: {len(pending_samples)} sample(s) "
                 f"({sum(1 for r in status_rows if r['overall']=='RESUME')} RESUME, "
                 f"{sum(1 for r in status_rows if r['overall']==_ST_TODO)} TODO)")
    # ──────────────────────────────────────────────────────────────

    # 2. Generate Snakefile (pending samples only, with priority map)
    smk_dir = out_root / ".snakemake_generated"
    smk_dir.mkdir(parents=True, exist_ok=True)
    snakefile_path = smk_dir / "Snakefile"
    snakefile_text = generate_snakefile(
        pending_samples, out_root, params, cmd, status_rows=status_rows)
    snakefile_path.write_text(snakefile_text, encoding="utf-8")
    logging.info(f"[SMK] Generated Snakefile: {snakefile_path}")
    logging.info(f"[SMK] Samples in Snakefile: {len(pending_samples)}  Subcommand: {cmd}")

    # ── v6.21: Always auto-unlock for embedded Snakemake ─────────
    # The whole point of wrapping Snakemake in Python is to spare users from
    # manual lock handling.  Therefore we *always* perform a safe unlock pass
    # before execution, even when no lock is currently visible.  This handles
    # stale lockfiles, partially-cleaned lock dirs, and bookkeeping drift.
    smk_unlock_flag = True
    lock_found = _smk_is_locked(out_root)

    if lock_found:
        logging.warning(
            "[SMK_UNLOCK] Stale Snakemake lock detected in "
            f"{out_root / '.snakemake' / 'locks'} "
            "— auto-unlocking before run.")
    else:
        logging.info(
            "[SMK_UNLOCK] Pre-run unlock check enabled by default "+
            "for embedded Snakemake backend.")
    rc_unlock = _smk_unlock(smk_exe, snakefile_path, out_root)
    if rc_unlock != 0:
        logging.error("[SMK_UNLOCK] Automatic unlock failed; aborting run.")
        return 1
    # ──────────────────────────────────────────────────────────────

    # 3. Determine conservative parallelism for embedded Snakemake
    user_max_parallel = max(1, int(params.get("max_parallel", 8)))
    eff_parallel = _auto_smk_spades_parallel(params, cmd)
    per_job_threads = max(
        int(params["threads_fastp"]),
        int(params["threads_spades"]),
        int(params["threads_busco"]),
    )
    sys_cpus = os.cpu_count() or 1
    # Guard: smk_cores must be at least 1
    smk_cores = max(1, min(eff_parallel * per_job_threads, sys_cpus))
    logging.info(
        f"[SMK] Parallel planning: user max_parallel={user_max_parallel}, "
        f"effective_parallel={eff_parallel}, per_job_threads={per_job_threads}, "
        f"cores={smk_cores}")

    # 4. Build snakemake command
    smk_cmd: List[str] = [
        smk_exe,
        "--snakefile", str(snakefile_path),
        "--cores", str(smk_cores),
        "--jobs", str(eff_parallel),
        "--rerun-incomplete",
        "--printshellcmds",
    ]
    if params.get("keep_going", True):
        smk_cmd.append("--keep-going")
    if dryrun:
        smk_cmd.append("--dryrun")

    # ── Auto --resources mem_mb to prevent OOM / server crash ────
    # Without this, Snakemake ignores resources.mem_mb and can launch
    # N SPAdes jobs simultaneously, each using 32GB → system crash.
    # Formula: total_mem_mb = system_ram * 85%
    # User can override with --smk_extra --resources mem_mb=XXXXX
    user_has_resources = any("--resources" in a for a in
                             list(params.get("smk_extra", []) or []))
    if not user_has_resources and not dryrun:
        try:
            with open("/proc/meminfo") as _f:
                for _ln in _f:
                    if _ln.startswith("MemTotal"):
                        _total_kb = int(_ln.split()[1])
                        _safe_mb = int(_total_kb / 1024 * 0.85)
                        smk_cmd += ["--resources", f"mem_mb={_safe_mb}"]
                        logging.info(
                            f"[SMK] Auto --resources mem_mb={_safe_mb} "
                            f"(85% of {_total_kb//1024} MB RAM)")
                        break
        except Exception:
            # fallback: 60GB safe limit
            smk_cmd += ["--resources", "mem_mb=61440"]
            logging.info("[SMK] Auto --resources mem_mb=61440 (fallback)")
    # ─────────────────────────────────────────────────────────────

    # ── v6.7: --prioritize for RESUME samples ────────────────────
    # Snakemake --prioritize bumps all jobs contributing to a target.
    # This ensures RESUME samples run before fresh TODO samples.
    resume_sids = [r["sid"] for r in status_rows if r["overall"] == "RESUME"]
    if resume_sids:
        _rf = cmd in ("all","fastp") and params.get("use_fastp",True)
        _rs = cmd in ("all","spades","spades2busco")
        _rb = cmd in ("all","busco","spades2busco")
        _sfx = params.get("suffix_spades_out","") or ""
        _pri_targets: List[str] = []
        for _sid in resume_sids[:50]:
            if _rb:
                _pri_targets.append(
                    str(out_root/params["busco_dir"]/_sid/".busco.done"))
            elif _rs:
                _pri_targets.append(
                    str(out_root/params["spades_dir"]/f"{_sid}{_sfx}"/".spades.done"))
            elif _rf:
                _pri_targets.append(
                    str(out_root/params["fastp_dir"]/_sid/".fastp.done"))
        if _pri_targets:
            smk_cmd.append("--prioritize")
            smk_cmd.extend(_pri_targets)
            logging.info(
                f"[SMK] --prioritize {len(_pri_targets)} RESUME samples")
    # ──────────────────────────────────────────────────────────────

    # ── v6.4: --smk_extra passthrough ────────────────────────────
    smk_extra_args: List[str] = list(params.get("smk_extra", []) or [])
    if smk_extra_args:
        logging.info(f"[SMK] Extra args: {' '.join(smk_extra_args)}")
        smk_cmd.extend(smk_extra_args)
    # ──────────────────────────────────────────────────────────────

    # 5. Execute snakemake (stream output)
    log_path = out_root / "snakemake.log"
    rc = _smk_run(smk_cmd, out_root, log_path)

    # 6. Handle lock error/race condition by auto-unlocking and retrying once.
    #    Even after a pre-run unlock, another process or an interrupted prior
    #    run can recreate the lock between check and execution.
    if rc != 0 and _smk_is_locked(out_root):
        logging.warning(
            "[SMK_UNLOCK] Lock appeared during/after run — auto-unlocking and retrying once.")
        if _smk_unlock(smk_exe, snakefile_path, out_root) == 0:
            rc = _smk_run(smk_cmd, out_root, log_path)

    if rc == 0:
        logging.info("[SMK] Snakemake finished successfully")
        if not dryrun:
            write_summary_tsv(samples, out_root, params)
    else:
        logging.error(
            f"[SMK] Snakemake exited with rc={rc}.\n"
            f"  Full log: {log_path}\n"
            f"  If error is 'Directory is locked', re-run with --smk_unlock")

    return 0 if rc == 0 else 1


# ==============================================================
# 12c. Background / daemonise (NEW v6)
# ==============================================================

def launch_background(out_root: Path) -> int:
    out_root.mkdir(parents=True, exist_ok=True)
    log_file = out_root / "pipeline.log"
    pid_file = out_root / "pipeline.pid"

    new_argv: List[str] = []
    for i, arg in enumerate(sys.argv):
        if arg in ("--background", "--bg", "--auto_bg"):
            continue
        new_argv.append(arg)

    # FIX v6.3: auto-inject --yes because background process has stdin=/dev/null
    # Without --yes, ask_user_confirmation() gets EOFError and silently aborts.
    if "--yes" not in new_argv and "-y" not in new_argv:
        new_argv.append("--yes")

    # Reconstruct command: always [sys.executable, script_path, ...args]
    # sys.argv[0] may be the script path or sys.executable depending on invocation.
    if new_argv and new_argv[0] == sys.executable:
        # e.g. python /path/to/script.py ... → already has interpreter
        cmd = new_argv
    else:
        # e.g. ./script.py ... or just script.py ...
        cmd = [sys.executable] + new_argv

    with open(str(log_file), "a", encoding="utf-8") as lf:
        lf.write(f"\n{'='*60}\n")
        lf.write(f"Background launch: {' '.join(cmd)}\n")
        lf.write(f"{'='*60}\n\n")
        lf.flush()
        proc = subprocess.Popen(
            cmd,
            stdout=lf,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            start_new_session=True,
        )

    pid_file.write_text(str(proc.pid) + "\n", encoding="utf-8")

    print(f"Pipeline launched in background.")
    print(f"  PID      : {proc.pid}")
    print(f"  Log      : {log_file}")
    print(f"  PID file : {pid_file}")
    print()
    print(f"Monitor progress:")
    print(f"  tail -f {log_file}")
    print(f"Stop:")
    print(f"  kill {proc.pid}")
    return 0


# ==============================================================
# 13. Clean subcommand
# ==============================================================

def _list_sample_dirs(root: Path) -> List[Path]:
    return sorted(p for p in root.iterdir() if p.is_dir()) if root.exists() else []


def clean_spades_batch(out_root: Path, spades_dir_name: str,
                       data_type: str, keep_patterns: List[str],
                       only_samples: Optional[List[str]],
                       force: bool, dryrun: bool) -> int:
    spades_root = out_root / spades_dir_name
    if not spades_root.exists():
        logging.error(f"[CLEAN] spades dir not found: {spades_root}")
        return 2
    targets = _list_sample_dirs(spades_root)
    if only_samples:
        allow = {safe_name(x) for x in only_samples}
        targets = [p for p in targets
                   if any(safe_name(p.name) == a
                          or safe_name(p.name).startswith(a) for a in allow)]
    n_clean = n_skip = 0
    for sdir in targets:
        ok_files = spades_success_files(sdir, data_type)
        if not force and not ok_files and not done_ok(sdir / ".spades.done"):
            n_skip += 1
            continue
        logging.info(f"[CLEAN_SPADES] {sdir.name}")
        if not dryrun:
            spades_clean_dir(sdir, keep_patterns)
        n_clean += 1
    logging.info(f"[CLEAN_SPADES_DONE] cleaned={n_clean}  skipped={n_skip}")
    return 0


def clean_fastp_batch(out_root: Path, fastp_dir_name: str,
                      only_samples: Optional[List[str]],
                      force: bool, dryrun: bool) -> int:
    root = out_root / fastp_dir_name
    if not root.exists():
        logging.error(f"[CLEAN] fastp dir not found: {root}")
        return 2
    targets = _list_sample_dirs(root)
    if only_samples:
        allow = {safe_name(x) for x in only_samples}
        targets = [p for p in targets if safe_name(p.name) in allow]
    n = 0
    for sdir in targets:
        if not force and not done_ok(sdir / ".fastp.done"):
            continue
        for pat in ["*.tmp", "*.temp", "tmp", ".tmp"]:
            for p in sdir.glob(pat):
                if dryrun:
                    continue
                if p.is_dir():
                    shutil.rmtree(p, ignore_errors=True)
                elif p.is_file():
                    try:
                        p.unlink()
                    except Exception:
                        pass
        n += 1
        logging.info(f"[CLEAN_FASTP] {sdir.name}")
    logging.info(f"[CLEAN_FASTP_DONE] cleaned={n}")
    return 0


def clean_busco_batch(out_root: Path, busco_dir_name: str,
                      only_samples: Optional[List[str]],
                      force: bool, dryrun: bool) -> int:
    root = out_root / busco_dir_name
    if not root.exists():
        logging.error(f"[CLEAN] busco dir not found: {root}")
        return 2
    targets = _list_sample_dirs(root)
    if only_samples:
        allow = {safe_name(x) for x in only_samples}
        targets = [p for p in targets if safe_name(p.name) in allow]
    n = 0
    for sdir in targets:
        if not force and not done_ok(sdir / ".busco.done"):
            continue
        for dname in ["tmp", "augustus_output", "metaeuk_output"]:
            d = sdir / dname
            if d.is_dir() and not dryrun:
                shutil.rmtree(d, ignore_errors=True)
        for d in list(sdir.glob("run_*/tmp")):
            if d.is_dir() and not dryrun:
                shutil.rmtree(d, ignore_errors=True)
        n += 1
        logging.info(f"[CLEAN_BUSCO] {sdir.name}")
    logging.info(f"[CLEAN_BUSCO_DONE] cleaned={n}")
    return 0


# ==============================================================
# 14. CLI  (v6.1: added --yes/-y)
# ==============================================================

_SUFFIX_HELP = textwrap.dedent("""\
    Fastq R1/R2 suffix specification. Accepts:
      auto                   - scan dir, detect automatically (DEFAULT)
      PRESET_NAME            - vendor preset: illumina, bgi, mgi, novogene,
                               macrogen, generic
      R1_SUF R2_SUF          - explicit pair, e.g.  _R1.fq.gz _R2.fq.gz
      PRESET1 PRESET2        - combine presets: illumina bgi
      R1a R2a R1b R2b ...    - multiple explicit pairs (even count)
    Examples:
      --suffix_fq auto
      --suffix_fq illumina
      --suffix_fq illumina bgi
      --suffix_fq _R1_001.fastq.gz _R2_001.fastq.gz
      --suffix_fq _1.fq.gz _2.fq.gz  _R1.fastq.gz _R2.fastq.gz""")


def _add_io(sp: argparse.ArgumentParser, input_required: bool = False) -> None:
    io = sp.add_argument_group("I/O")
    grp = io.add_mutually_exclusive_group(required=input_required)
    grp.add_argument("-i", "--input_dir", default=None,
                     help="raw reads directory\n  Example: -i /data/reads")
    grp.add_argument("--sample_table", default=None, metavar="TSV",
                     help="sample table TSV\n  Example: --sample_table samples.tsv")
    io.add_argument("-o", "--out_root", default="pipeline_out",
                    help="output root [%(default)s]")
    io.add_argument("--fastp_dir", default="01.fastp")
    io.add_argument("--spades_dir", default="02.spades")
    io.add_argument("--busco_dir", default="03.busco")
    io.add_argument("--summary_tsv_name", default="summary_report.tsv")


def _add_scan(sp: argparse.ArgumentParser) -> None:
    g = sp.add_argument_group("Reads scanning")
    g.add_argument("--fq_position", choices=["subdir", "flat"], default="flat",
                   help="flat: all fastq in --input_dir; subdir: per-sample dirs [%(default)s]")
    g.add_argument("--suffix_fq", nargs="*", default=["auto"], metavar="SPEC",
                   help=_SUFFIX_HELP)
    g.add_argument("--subdir_list", nargs="*", default=None,
                   help="restrict to these subdirs (space/comma/file)")


def _add_mode(sp: argparse.ArgumentParser) -> None:
    g = sp.add_argument_group("Data type / assembly mode")
    g.add_argument("--data_type", choices=["wgs", "rna"], default="wgs")
    g.add_argument("--assembly_mode", default="",
                   help="SPAdes mode: " + ", ".join(sorted(VALID_ASM_MODES)))
    g.add_argument("--which_fasta_to_busco", "-wfb", default="AUTO")


def _add_runner(sp: argparse.ArgumentParser) -> None:
    g = sp.add_argument_group("Runner")
    g.add_argument("--dryrun", action="store_true")
    g.add_argument("--max_parallel", type=int, default=8)
    g.add_argument("--keep_going", action=argparse.BooleanOptionalAction, default=True)
    g.add_argument("--skip_samples", nargs="*", default=None, metavar="ID_OR_FILE")
    g.add_argument("--no_preflight", action="store_true")
    g.add_argument("--resource_cap_pct", type=float, default=90.0, metavar="PCT")
    # --- v6 additions ---
    g.add_argument("--smk", "--snakemake", dest="use_snakemake",
                   action="store_true", default=False,
                   help="use snakemake backend (requires snakemake in PATH)")
    g.add_argument("--smk_extra", nargs=argparse.REMAINDER, default=[],
                   metavar="ARG",
                   help="extra args forwarded to snakemake (must be last)")
    g.add_argument("--smk_unlock", "--unlock", dest="smk_unlock",
                   action="store_true", default=False,
                   help="optional legacy flag; auto-unlock is already built in by default")
    g.add_argument("--smk_mem_mb_fastp",  type=int, default=8000,  metavar="MB",
                   help="mem_mb for fastp rule [%(default)s]")
    g.add_argument("--smk_mem_mb_spades", type=int, default=0,     metavar="MB",
                   help="mem_mb for SPAdes rule [auto=spades_mem_gb*1024]")
    g.add_argument("--smk_mem_mb_busco",  type=int, default=16000, metavar="MB",
                   help="mem_mb for BUSCO rule [%(default)s]")
    g.add_argument("--smk_disk_mb",       type=int, default=200000, metavar="MB",
                   help="disk_mb for all rules [%(default)s]")
    g.add_argument("--background", "--bg", dest="background",
                   action="store_true", default=False,
                   help="run pipeline in background (immune to SSH disconnect)")
    # --- v6.1 addition ---
    g.add_argument("--yes", "-y", dest="auto_yes",
                   action="store_true", default=False,
                   help="skip interactive confirmation, run directly")
    # --- v6.3 additions ---
    g.add_argument("--auto_bg", dest="auto_bg",
                   action="store_true", default=False,
                   help="run first --pilot_n samples, then re-launch rest in background")
    g.add_argument("--pilot_n", type=int, default=5, metavar="N",
                   help="pilot samples to run before auto-bg [%(default)s]")


def _add_tools(sp: argparse.ArgumentParser) -> None:
    g = sp.add_argument_group("Tool execution")
    g.add_argument("--conda_exe", default="conda",
                   help="conda/mamba executable [%(default)s]")
    # FIX v6.4: default=None so preflight does NOT require conda unless
    # user explicitly sets --conda_env_*.  Old defaults ("fastp","spades","busco")
    # forced a conda check on every run even when tools are on PATH directly.
    g.add_argument("--conda_env_fastp",  default=None, metavar="ENV",
                   help="conda env containing fastp (omit to use PATH)")
    g.add_argument("--conda_env_spades", default=None, metavar="ENV",
                   help="conda env containing spades.py (omit to use PATH)")
    g.add_argument("--conda_env_busco",  default=None, metavar="ENV",
                   help="conda env containing busco (omit to use PATH)")
    g.add_argument("--fastp_exe",  default="fastp",     metavar="EXE")
    g.add_argument("--spades_exe", default="spades.py", metavar="EXE")
    g.add_argument("--busco_exe",  default="busco",     metavar="EXE")
    g.add_argument("--fastp_path",  default=None, metavar="PATH")
    g.add_argument("--spades_path", default=None, metavar="PATH")
    g.add_argument("--busco_path",  default=None, metavar="PATH")


def _add_fastp_opts(sp: argparse.ArgumentParser) -> None:
    g = sp.add_argument_group("fastp options")
    g.add_argument("--use_fastp", action=argparse.BooleanOptionalAction, default=True)
    g.add_argument("--threads_fastp", type=int, default=8)
    g.add_argument("--fastp_detect_adapter_for_pe",
                   action=argparse.BooleanOptionalAction, default=True)
    g.add_argument("--fastp_trim_poly_x",
                   action=argparse.BooleanOptionalAction, default=True)
    g.add_argument("--fastp_reads_to_process", type=int, default=0)
    g.add_argument("--fastp_extra", action="append", default=[])


def _add_spades_opts(sp: argparse.ArgumentParser) -> None:
    g = sp.add_argument_group("SPAdes options")
    g.add_argument("--threads_spades", type=int, default=16)
    g.add_argument("--spades_mem_gb", type=int, default=32)
    g.add_argument("--spades_extra", action="append", default=[])
    g.add_argument("--spades_clean", action=argparse.BooleanOptionalAction, default=True)
    g.add_argument("--keep_files", action="append", default=None)
    g.add_argument("--suffix_spades_out", default="", metavar="SUFFIX",
                   help="e.g. --suffix_spades_out _k77")


def _add_busco_opts(sp: argparse.ArgumentParser,
                    lineage_required: bool = True) -> None:
    g = sp.add_argument_group("BUSCO options")
    g.add_argument("--busco_lineage", required=lineage_required,
                   default="" if not lineage_required else None,
                   help="e.g. bacteria_odb10 OR /path/to/bacteria_odb10")
    g.add_argument("--threads_busco", type=int, default=12)
    g.add_argument("--busco_mode", default="")
    g.add_argument("--busco_offline", action=argparse.BooleanOptionalAction, default=False)
    g.add_argument("--busco_force", action=argparse.BooleanOptionalAction, default=False)
    g.add_argument("--busco_extra", action="append", default=[])
    g.add_argument("--busco_fasta_dir", default=None, metavar="DIR")


def _add_hidden_busco(sp: argparse.ArgumentParser) -> None:
    for flag, default in [("--busco_lineage", ""), ("--busco_mode", ""),
                          ("--busco_fasta_dir", None)]:
        sp.add_argument(flag, default=default, help=argparse.SUPPRESS)
    for flag in ["--busco_offline", "--busco_force"]:
        sp.add_argument(flag, action=argparse.BooleanOptionalAction,
                        default=False, help=argparse.SUPPRESS)
    sp.add_argument("--busco_extra", action="append", default=[], help=argparse.SUPPRESS)
    sp.add_argument("--threads_busco", type=int, default=12, help=argparse.SUPPRESS)


def _add_hidden_fastp_spades(sp: argparse.ArgumentParser) -> None:
    for flag, default in [("--use_fastp", False),
                          ("--fastp_detect_adapter_for_pe", True),
                          ("--fastp_trim_poly_x", True),
                          ("--spades_clean", True)]:
        sp.add_argument(flag, action=argparse.BooleanOptionalAction,
                        default=default, help=argparse.SUPPRESS)
    for flag, dv in [("--threads_fastp", 8), ("--fastp_reads_to_process", 0),
                     ("--threads_spades", 16), ("--spades_mem_gb", 32)]:
        sp.add_argument(flag, type=int, default=dv, help=argparse.SUPPRESS)
    for flag in ["--fastp_extra", "--spades_extra"]:
        sp.add_argument(flag, action="append", default=[], help=argparse.SUPPRESS)
    sp.add_argument("--keep_files", action="append", default=None, help=argparse.SUPPRESS)
    sp.add_argument("--suffix_spades_out", default="", help=argparse.SUPPRESS)


def build_parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(
        prog="fastp_spades_busco_v6_25.py",
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent("""\
        --- Quick-start examples ---

        Full pipeline (preview plan, then confirm):
          %(prog)s all -i /data/reads -o /results \\
              --suffix_fq illumina --busco_lineage bacteria_odb10

        Full pipeline (skip confirmation):
          %(prog)s all -i /data/reads -o /results \\
              --suffix_fq auto --busco_lineage bacteria_odb10 --yes

        Full pipeline (background, snakemake):
          %(prog)s all -i /data/reads -o /results \\
              --busco_lineage bacteria_odb10 --smk --bg --yes

        Full pipeline (local BUSCO lineage, auto-offline):
          %(prog)s all -i /data/reads -o /results \\
              --busco_lineage /db/busco/bacteria_odb10 --smk --yes

        Snakemake + cluster (SLURM):
          %(prog)s all -i /data/reads -o /results \\
              --busco_lineage metazoa_odb10 --smk --yes \\
              --smk_extra --cluster "sbatch -c {threads} --mem=32G -t 24:00:00" \\
                          --latency-wait 120 --jobs 50

        Unlock stale Snakemake lock (after crash/kill):
          %(prog)s all -i /data/reads -o /results \\
              --busco_lineage metazoa_odb10 --smk --yes --smk_unlock (legacy; auto-enabled)

        SPAdes only:
          %(prog)s spades -i /data/reads -o /results/assembly --yes

        BUSCO on existing assemblies:
          %(prog)s busco -o /results --busco_fasta_dir /assemblies \\
              --busco_lineage bacteria_odb10 --yes

        Check environment:
          %(prog)s check --conda_env_spades spades_env

        Show task status (no run):
          %(prog)s status -o /results
          %(prog)s status -o /results --show_done
          %(prog)s status -i /data/reads -o /results --suffix_fq illumina

        Snakemake with HPC resources:
          %(prog)s all -i /data/reads -o /results --busco_lineage metazoa_odb10 \\
              --smk --yes \\
              --smk_mem_mb_spades 65536 --smk_disk_mb 500000 \\
              --smk_extra --cluster "sbatch -c {threads} --mem={resources.mem_mb}M" \\
                          --latency-wait 120 --jobs 50

        --- Supported --suffix_fq presets ---
        """ + "\n".join(
            f"  {name:<12s}  {pairs[0][0]} / {pairs[0][1]}"
            + (f"  (+ {len(pairs)-1} more)" if len(pairs) > 1 else "")
            for name, pairs in sorted(_SUFFIX_PRESETS.items())
        )),
    )
    sub = root.add_subparsers(dest="cmd", metavar="SUBCOMMAND")

    # -- all --------------------------------------------------------
    sp = sub.add_parser("all", help="fastp -> SPAdes -> BUSCO (full)")
    _add_io(sp, True); _add_scan(sp); _add_mode(sp); _add_runner(sp)
    _add_tools(sp); _add_fastp_opts(sp); _add_spades_opts(sp)
    _add_busco_opts(sp, True)

    # -- fastp ------------------------------------------------------
    sp = sub.add_parser("fastp", help="run fastp only")
    _add_io(sp, True); _add_scan(sp); _add_runner(sp); _add_tools(sp)
    _add_fastp_opts(sp)
    sp.add_argument("--data_type", choices=["wgs", "rna"], default="wgs")
    _add_hidden_busco(sp)
    for flag, dv in [("--threads_spades", 16), ("--spades_mem_gb", 32)]:
        sp.add_argument(flag, type=int, default=dv, help=argparse.SUPPRESS)
    sp.add_argument("--spades_extra", action="append", default=[], help=argparse.SUPPRESS)
    sp.add_argument("--spades_clean", action=argparse.BooleanOptionalAction,
                    default=True, help=argparse.SUPPRESS)
    sp.add_argument("--keep_files", action="append", default=None, help=argparse.SUPPRESS)
    sp.add_argument("--assembly_mode", default="", help=argparse.SUPPRESS)
    sp.add_argument("--which_fasta_to_busco", default="AUTO", help=argparse.SUPPRESS)
    sp.add_argument("--suffix_spades_out", default="", help=argparse.SUPPRESS)

    # -- spades -----------------------------------------------------
    sp = sub.add_parser("spades", help="run SPAdes only (output into -o/<sample>/)")
    _add_io(sp, True); _add_scan(sp); _add_mode(sp); _add_runner(sp)
    _add_tools(sp); _add_spades_opts(sp); _add_hidden_busco(sp)
    for flag, default in [("--fastp_detect_adapter_for_pe", True),
                          ("--fastp_trim_poly_x", True)]:
        sp.add_argument(flag, action=argparse.BooleanOptionalAction,
                        default=default, help=argparse.SUPPRESS)
    sp.add_argument("--use_fastp", action=argparse.BooleanOptionalAction,
                    default=False, help=argparse.SUPPRESS)
    for flag, dv in [("--threads_fastp", 8), ("--fastp_reads_to_process", 0)]:
        sp.add_argument(flag, type=int, default=dv, help=argparse.SUPPRESS)
    sp.add_argument("--fastp_extra", action="append", default=[], help=argparse.SUPPRESS)

    # -- spades2busco -----------------------------------------------
    sp = sub.add_parser("spades2busco", help="SPAdes -> BUSCO (no fastp)")
    _add_io(sp, True); _add_scan(sp); _add_mode(sp); _add_runner(sp)
    _add_tools(sp); _add_spades_opts(sp); _add_busco_opts(sp, True)
    for flag, default in [("--fastp_detect_adapter_for_pe", True),
                          ("--fastp_trim_poly_x", True)]:
        sp.add_argument(flag, action=argparse.BooleanOptionalAction,
                        default=default, help=argparse.SUPPRESS)
    sp.add_argument("--use_fastp", action=argparse.BooleanOptionalAction,
                    default=False, help=argparse.SUPPRESS)
    for flag, dv in [("--threads_fastp", 8), ("--fastp_reads_to_process", 0)]:
        sp.add_argument(flag, type=int, default=dv, help=argparse.SUPPRESS)
    sp.add_argument("--fastp_extra", action="append", default=[], help=argparse.SUPPRESS)

    # -- busco ------------------------------------------------------
    sp = sub.add_parser("busco", help="run BUSCO only (output into -o/<sample>/)")
    _add_io(sp, False); _add_scan(sp); _add_mode(sp); _add_runner(sp)
    _add_tools(sp); _add_busco_opts(sp, True); _add_hidden_fastp_spades(sp)

    # -- check ------------------------------------------------------
    sp = sub.add_parser("check", help="verify tool availability")
    _add_tools(sp)
    sp.add_argument("-o", "--out_root", default="pipeline_out", help=argparse.SUPPRESS)
    sp.add_argument("--stages", nargs="*", choices=["fastp", "spades", "busco"],
                    default=["fastp", "spades", "busco"])
    sp.add_argument("--use_fastp", action=argparse.BooleanOptionalAction,
                    default=True, help=argparse.SUPPRESS)
    sp.add_argument("--resource_cap_pct", type=float, default=90.0, help=argparse.SUPPRESS)

    # -- status -----------------------------------------------------
    sp = sub.add_parser("status",
                        help="scan existing outputs and show completion table (no run)")
    _add_io(sp, False)
    _add_scan(sp)
    _add_mode(sp)
    _add_tools(sp)
    _add_fastp_opts(sp)
    _add_spades_opts(sp)
    _add_busco_opts(sp, lineage_required=False)
    sp.add_argument("--max_rows", type=int, default=200,
                    help="max rows to print (use 0 for unlimited) [%(default)s]")
    sp.add_argument("--show_done", action="store_true", default=False,
                    help="include already-complete samples in the table")
    sp.add_argument(
        "--pipeline_cmd",
        choices=["auto", "all", "fastp", "spades", "busco", "spades2busco"],
        default="auto",
        help=(
            "Which pipeline stages to show in the status table.  "
            "'auto' (default) infers from existing output directories: "
            "presence of 01.fastp/ -> fastp; 02.spades/ -> spades; "
            "03.busco/ -> busco.  Override with 'all', 'spades', 'busco', "
            "'spades2busco', or 'fastp' to force a specific view."
        ),
    )
    # hidden stubs so make_params() doesn't fail
    sp.add_argument("--resource_cap_pct", type=float, default=90.0,
                    help=argparse.SUPPRESS)
    sp.add_argument("--max_parallel", type=int, default=8, help=argparse.SUPPRESS)
    sp.add_argument("--keep_going",   action=argparse.BooleanOptionalAction,
                    default=True, help=argparse.SUPPRESS)
    sp.add_argument("--skip_samples", nargs="*", default=None,
                    help=argparse.SUPPRESS)
    sp.add_argument("--no_preflight", action="store_true", help=argparse.SUPPRESS)
    sp.add_argument("--dryrun",       action="store_true", help=argparse.SUPPRESS)
    sp.add_argument("--smk", "--snakemake", dest="use_snakemake",
                    action="store_true", default=False, help=argparse.SUPPRESS)
    sp.add_argument("--background", "--bg", dest="background",
                    action="store_true", default=False, help=argparse.SUPPRESS)
    sp.add_argument("--yes", "-y", dest="auto_yes",
                    action="store_true", default=False, help=argparse.SUPPRESS)
    sp.add_argument("--auto_bg",  action="store_true", default=False,
                    help=argparse.SUPPRESS)
    sp.add_argument("--pilot_n",  type=int, default=5, help=argparse.SUPPRESS)
    sp.add_argument("--smk_extra", nargs=argparse.REMAINDER, default=[],
                    help=argparse.SUPPRESS)
    sp.add_argument("--smk_unlock", action="store_true", default=False,
                    help=argparse.SUPPRESS)
    sp.add_argument("--smk_mem_mb_fastp",  type=int, default=8000,
                    help=argparse.SUPPRESS)
    sp.add_argument("--smk_mem_mb_spades", type=int, default=0,
                    help=argparse.SUPPRESS)
    sp.add_argument("--smk_mem_mb_busco",  type=int, default=16000,
                    help=argparse.SUPPRESS)
    sp.add_argument("--smk_disk_mb",       type=int, default=200000,
                    help=argparse.SUPPRESS)

    # -- clean ------------------------------------------------------
    sp = sub.add_parser("clean", help="remove intermediate files")
    sp.add_argument("-o", "--out_root", required=True)
    sp.add_argument("--fastp_dir", default="01.fastp")
    sp.add_argument("--spades_dir", default="02.spades")
    sp.add_argument("--busco_dir", default="03.busco")
    sp.add_argument("--data_type", choices=["wgs", "rna"], default="wgs")
    sp.add_argument("--keep_files", action="append", default=None)
    sp.add_argument("--clean_stage", choices=["spades", "fastp", "busco", "all"],
                    default="spades")
    sp.add_argument("--clean_samples", nargs="*", default=None)
    sp.add_argument("--clean_force", action="store_true")
    sp.add_argument("--dryrun", action="store_true")
    sp.add_argument("--resource_cap_pct", type=float, default=90.0, help=argparse.SUPPRESS)

    return root


def make_params(args) -> dict:
    G = getattr
    user_keep = parse_keep_files(G(args, "keep_files", None))
    keep_pats = user_keep if user_keep else default_keep_patterns()
    for x in [".spades.done", ".spades.fail", "*.log", "*.txt", "*.yaml", "*.yml"]:
        if x not in keep_pats:
            keep_pats.append(x)
    return {
        "conda_exe":       which_or_none(G(args, "conda_exe", "conda")) or G(args, "conda_exe", "conda"),
        "conda_env_fastp": G(args, "conda_env_fastp", None),
        "conda_env_spades": G(args, "conda_env_spades", None),
        "conda_env_busco": G(args, "conda_env_busco", None),
        "fastp_exe":       G(args, "fastp_exe", "fastp"),
        "spades_exe":      G(args, "spades_exe", "spades.py"),
        "busco_exe":       G(args, "busco_exe", "busco"),
        "fastp_path":      G(args, "fastp_path", None),
        "spades_path":     G(args, "spades_path", None),
        "busco_path":      G(args, "busco_path", None),
        "fastp_dir":       G(args, "fastp_dir", "01.fastp"),
        "spades_dir":      G(args, "spades_dir", "02.spades"),
        "busco_dir":       G(args, "busco_dir", "03.busco"),
        "summary_tsv_name": G(args, "summary_tsv_name", "summary_report.tsv"),
        "data_type":       G(args, "data_type", "wgs"),
        "assembly_mode":   normalize_assembly_mode(G(args, "assembly_mode", "")),
        "which_fasta_to_busco": G(args, "which_fasta_to_busco", "AUTO"),
        "use_fastp":       bool(G(args, "use_fastp", True)),
        "threads_fastp":   int(G(args, "threads_fastp", 8)),
        "fastp_reads_to_process": int(G(args, "fastp_reads_to_process", 0) or 0),
        "fastp_detect_adapter_for_pe": bool(G(args, "fastp_detect_adapter_for_pe", True)),
        "fastp_trim_poly_x": bool(G(args, "fastp_trim_poly_x", True)),
        "fastp_extra":     list(G(args, "fastp_extra", []) or []),
        "threads_spades":  int(G(args, "threads_spades", 16)),
        "spades_mem_gb":   int(G(args, "spades_mem_gb", 32)),
        "spades_extra":    list(G(args, "spades_extra", []) or []),
        "spades_clean":    bool(G(args, "spades_clean", True)),
        "keep_files":      keep_pats,
        "suffix_spades_out": (G(args, "suffix_spades_out", "") or ""),
        "direct_spades_output": False,
        "direct_busco_output": False,
        "busco_lineage":   (G(args, "busco_lineage", "") or ""),
        "threads_busco":   int(G(args, "threads_busco", 12)),
        "busco_mode":      (G(args, "busco_mode", "") or "").strip(),
        "busco_offline":   bool(G(args, "busco_offline", False)),
        "busco_force":     bool(G(args, "busco_force", False)),
        "busco_extra":     list(G(args, "busco_extra", []) or []),
        "busco_fasta_dir": G(args, "busco_fasta_dir", None),
        "keep_going":      bool(G(args, "keep_going", True)),
        "max_parallel":    int(G(args, "max_parallel", 8)),
        "skip_samples_raw": G(args, "skip_samples", None),
        "resource_cap_pct": float(G(args, "resource_cap_pct", 90.0)),
        # v6.4: snakemake passthrough
        "smk_extra":  list(G(args, "smk_extra", []) or []),
        "smk_unlock": bool(G(args, "smk_unlock", False)),
        "auto_bg":    bool(G(args, "auto_bg",    False)),
        "pilot_n":    int(G(args,  "pilot_n",    5)),
        # v6.5: snakemake resource declarations
        "smk_mem_mb_fastp":  max(1, int(G(args, "smk_mem_mb_fastp",  8000)  or 8000)),
        "smk_mem_mb_spades": max(1, int(G(args, "smk_mem_mb_spades", 0) or (int(G(args, "spades_mem_gb", 32) or 32) * 1024))),
        "smk_mem_mb_busco":  max(1, int(G(args, "smk_mem_mb_busco",  16000) or 16000)),
        "smk_disk_mb":       int(G(args, "smk_disk_mb",       200000) or 200000),
    }


# ==============================================================
# 15. main
# ==============================================================

def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if not args.cmd:
        parser.print_help()
        return 0

    setup_root_logger()

    # ============================================================
    # v6: --background  (must come BEFORE any heavy work)
    # ============================================================
    if getattr(args, "background", False):
        out_root = Path(getattr(args, "out_root", "pipeline_out")).resolve()
        return launch_background(out_root)

    # == clean ==
    if args.cmd == "clean":
        out_root = Path(args.out_root).resolve()
        out_root.mkdir(parents=True, exist_ok=True)
        setup_root_logger(out_root / "pipeline.log")
        keep_pats = parse_keep_files(args.keep_files) or default_keep_patterns()
        # Mirror the safety patch from make_params(): always preserve pipeline
        # stamps and log/txt files so Snakemake can find them after a clean.
        for _x in [".spades.done", ".spades.fail", "*.log", "*.txt", "*.yaml", "*.yml"]:
            if _x not in keep_pats:
                keep_pats.append(_x)
        only = (resolve_list_arg(args.clean_samples, out_root)
                if args.clean_samples else None)
        rc = 0
        if args.clean_stage in ("spades", "all"):
            rc = clean_spades_batch(out_root, args.spades_dir, args.data_type,
                                    keep_pats, only, args.clean_force, args.dryrun)
        if rc == 0 and args.clean_stage in ("fastp", "all"):
            rc = clean_fastp_batch(out_root, args.fastp_dir,
                                   only, args.clean_force, args.dryrun)
        if rc == 0 and args.clean_stage in ("busco", "all"):
            rc = clean_busco_batch(out_root, args.busco_dir,
                                   only, args.clean_force, args.dryrun)
        logging.info("[CLEAN_DONE]" if rc == 0 else f"[CLEAN_ERROR] rc={rc}")
        return rc

    # == check ==
    if args.cmd == "check":
        params = make_params(args)
        out_root = Path(args.out_root).resolve()
        stages: Set[str] = set(
            getattr(args, "stages", None) or ["fastp", "spades", "busco"])
        if not stages:
            stages = {"fastp", "spades", "busco"}
        use_smk = bool(getattr(args, "use_snakemake", False))
        return preflight_tools_only_check(stages, params, out_root, use_smk)

    # == status ==
    if args.cmd == "status":
        out_root = Path(args.out_root).resolve()
        setup_root_logger()  # console only for status
        params = make_params(args)
        params["_out_root"] = str(out_root)
        max_rows  = int(getattr(args, "max_rows",  200) or 200)
        show_done = bool(getattr(args, "show_done", False))

        # Collect samples (best-effort, require_fastq=False for busco-only dirs)
        raw_suffix = getattr(args, "suffix_fq", ["auto"]) or ["auto"]
        try:
            suffix_pairs = resolve_suffix_fq(raw_suffix)
        except ValueError as e:
            logging.error(f"[ARG_ERROR] {e}")
            return 1

        # Determine which pipeline stages to use for the status table.
        # BUG-FIX 17: was hardcoded to "all", which treated every stage as
        # active and dragged "overall" to TODO when stages like fastp were
        # never run.  Now auto-detected from directory structure after
        # sample collection; can be overridden with --pipeline_cmd.
        _explicit_pipeline_cmd = getattr(args, "pipeline_cmd", "auto") or "auto"
        if getattr(args, "input_dir", None) or getattr(args, "sample_table", None):
            try:
                samples, _ = collect_samples(args, suffix_pairs, require_fastq=False)
            except RuntimeError as e:
                logging.error(str(e)); return 1
        else:
            # Auto-scan output directories for sample IDs.
            # Strategy (in priority order):
            #   1. Nested SPAdes dir  (standard full-pipeline: out_root/spades_dir/<sid>/)
            #   2. Nested BUSCO dir   (busco-only: out_root/busco_dir/<sid>/)
            #   3. Direct out_root    (spades/busco subcommand: out_root/<sid>/)
            # In every case: strip --suffix_spades_out and skip pipeline/hidden dirs.
            suffix = params.get("suffix_spades_out", "") or ""
            pipeline_dir_names: Set[str] = {
                params.get("fastp_dir",  "01.fastp"),
                params.get("busco_dir",  "03.busco"),
                params.get("spades_dir", "02.spades"),
                ".snakemake", ".snakemake_generated", "__pycache__",
            }
            seen_sids: Set[str] = set()
            samples = []

            def _scan_dir_for_samples(dirpath: Path) -> None:
                """Append Sample stubs from *dirpath*, deduplicating by sid."""
                if not dirpath.is_dir():
                    return
                for d in sorted(dirpath.iterdir()):
                    if not d.is_dir():
                        continue
                    dname = d.name
                    # Skip hidden dirs and known pipeline infrastructure dirs
                    if dname.startswith(".") or dname in pipeline_dir_names:
                        continue
                    sid = (dname[: -len(suffix)]
                           if suffix and dname.endswith(suffix) else dname)
                    if sid and sid not in seen_sids:
                        seen_sids.add(sid)
                        samples.append(Sample(sid=sid, r1="", r2=""))

            # 1. Nested SPAdes output directory
            _scan_dir_for_samples(out_root / params.get("spades_dir", "02.spades"))
            # 2. Nested BUSCO output directory
            _scan_dir_for_samples(out_root / params.get("busco_dir", "03.busco"))
            # 3. Fall back to direct out_root scan (spades/busco subcommand mode)
            if not samples:
                _scan_dir_for_samples(out_root)

            if not samples:
                logging.error(
                    f"[STATUS] No samples found.\n"
                    f"  Tried: {out_root / params.get('spades_dir','02.spades')}\n"
                    f"         {out_root / params.get('busco_dir','03.busco')}\n"
                    f"         {out_root}  (direct output fallback)\n"
                    f"  Provide -i / --sample_table, or run from the pipeline output dir.")
                return 1

        # ── Derive status_cmd ──────────────────────────────────────────
        # If the user supplied an explicit --pipeline_cmd, honour it.
        # Otherwise infer from which stage output directories are present so
        # that stage columns not relevant to this particular run are hidden and
        # "overall" is computed correctly.
        if _explicit_pipeline_cmd != "auto":
            status_cmd = _explicit_pipeline_cmd
        else:
            _fastp_d  = out_root / params.get("fastp_dir",  "01.fastp")
            _spades_d = out_root / params.get("spades_dir", "02.spades")
            _busco_d  = out_root / params.get("busco_dir",  "03.busco")

            def _dir_has_samples(d: Path) -> bool:
                """True if *d* exists and contains at least one non-hidden subdir."""
                if not d.is_dir():
                    return False
                try:
                    return any(
                        child.is_dir() and not child.name.startswith(".")
                        for child in d.iterdir()
                    )
                except PermissionError:
                    return False

            _has_fastp_dir  = _dir_has_samples(_fastp_d)
            _has_spades_dir = _dir_has_samples(_spades_d)
            _has_busco_dir  = _dir_has_samples(_busco_d)

            # For direct-output mode (spades/busco subcommand), check stamp
            # files under the first discovered sample directory.
            _has_direct_spades = _has_direct_busco = False
            if not _has_spades_dir and not _has_busco_dir and samples:
                _first_sid = samples[0].sid
                _has_direct_spades = done_ok(out_root / _first_sid / ".spades.done")
                _has_direct_busco  = done_ok(out_root / _first_sid / ".busco.done")

            _has_spades = _has_spades_dir or _has_direct_spades
            _has_busco  = _has_busco_dir  or _has_direct_busco
            _use_fp     = params.get("use_fastp", True)

            if _has_fastp_dir and _use_fp and _has_spades and _has_busco:
                status_cmd = "all"
            elif _has_spades and _has_busco and not (_has_fastp_dir and _use_fp):
                status_cmd = "spades2busco"
            elif _has_spades and not _has_busco:
                status_cmd = "spades"
            elif _has_busco and not _has_spades:
                status_cmd = "busco"
            elif _has_fastp_dir and _use_fp:
                status_cmd = "fastp"
            else:
                status_cmd = "all"   # nothing detected yet → show everything

        logging.info(
            f"[STATUS] pipeline_cmd={status_cmd!r}  "
            f"(source={'--pipeline_cmd' if _explicit_pipeline_cmd != 'auto' else 'auto-detected'})"
        )

        rows = scan_all_status(samples, out_root, params, status_cmd)

        if not show_done:
            display_rows = [r for r in rows if r["overall"] != _ST_DONE]
            n_hidden = len(rows) - len(display_rows)
        else:
            display_rows = rows
            n_hidden = 0

        if max_rows and max_rows > 0:
            _max = max_rows
        else:
            _max = len(display_rows)

        print_status_table(display_rows[:_max], status_cmd, params,
                           out_root, max_rows=_max)

        n_done    = sum(1 for r in rows if r["overall"] == _ST_DONE)
        n_resume  = sum(1 for r in rows if r["overall"] == "RESUME")
        n_todo    = sum(1 for r in rows if r["overall"] == _ST_TODO)
        print(f"  Full count: total={len(rows)}  done={n_done}  "
              f"partial={n_resume}  todo={n_todo}")
        if n_hidden:
            print(f"  (Hiding {n_hidden} DONE samples; use --show_done to include them)")

        # Per-stage counts
        for stage in ("fastp", "spades", "busco"):
            stage_done    = sum(1 for r in rows if r[stage] == _ST_DONE)
            stage_partial = sum(1 for r in rows if r[stage] == _ST_PARTIAL)
            stage_todo    = sum(1 for r in rows if r[stage] == _ST_TODO)
            stage_skip    = sum(1 for r in rows if r[stage] == _ST_SKIP)
            if stage_skip == len(rows):
                continue  # stage not in this pipeline run
            print(f"  {stage:<8s}: done={stage_done}  partial={stage_partial}"
                  f"  todo={stage_todo}")
        return 0

    # == pipeline subcommands ==
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "pipeline.log")
    logging.info(f"[START] cmd={args.cmd}  out_root={out_root}")

    params = make_params(args)
    params["_out_root"] = str(out_root)

    # ============================================================
    # v6: auto-detect local BUSCO lineage -> --offline
    # ============================================================
    lineage = params.get("busco_lineage", "")
    if lineage and is_local_busco_lineage(lineage):
        # v6.16: resolve to absolute path — Snakemake runs from out_root as CWD,
        # so relative paths like '../../data/lineage' would resolve incorrectly.
        lineage_p = Path(lineage)
        if not lineage_p.is_absolute():
            lineage_p = lineage_p.resolve()
            params["busco_lineage"] = str(lineage_p)
            logging.info(f"[CONFIG] BUSCO lineage resolved to: {lineage_p}")
        if not params["busco_offline"]:
            params["busco_offline"] = True
            logging.info(
                f"[CONFIG] BUSCO lineage is a local path: {params['busco_lineage']}\n"
                f"         -> auto-enabling --offline (no download)")

    # v6.16: resolve relative busco_fasta_dir — same CWD issue as lineage
    bfd = params.get("busco_fasta_dir")
    if bfd:
        bfd_p = Path(bfd)
        if not bfd_p.is_absolute():
            bfd_p = bfd_p.resolve()
            params["busco_fasta_dir"] = str(bfd_p)
            logging.info(f"[CONFIG] busco_fasta_dir resolved to: {bfd_p}")

    # -- resolve suffix_fq --
    raw_suffix = getattr(args, "suffix_fq", ["auto"]) or ["auto"]
    try:
        suffix_pairs = resolve_suffix_fq(raw_suffix)
    except ValueError as e:
        logging.error(f"[ARG_ERROR] {e}")
        return 1
    if suffix_pairs:
        for r1s, r2s in suffix_pairs[:3]:
            logging.info(f"[CONFIG] suffix_pair: {r1s}  /  {r2s}")
        if len(suffix_pairs) > 3:
            logging.info(f"[CONFIG] ... and {len(suffix_pairs) - 3} more patterns")

    # -- direct-output flags --
    if args.cmd == "spades":
        params["direct_spades_output"] = True
        logging.info("[CONFIG] direct_spades_output=True")
    if args.cmd == "busco":
        params["direct_busco_output"] = True
        logging.info("[CONFIG] direct_busco_output=True")

    # -- resource cap --
    cap_pct = params["resource_cap_pct"] / 100.0
    params = cap_resources(params, max_pct=cap_pct)

    # v6.3: full preflight runs after sample collection (below)
    #        so sample-level checks have all the data they need.

    if params.get("suffix_spades_out"):
        logging.info(f"[CONFIG] suffix_spades_out={params['suffix_spades_out']!r}")

    # == Collect samples ==
    samples: List[Sample] = []

    if args.cmd == "busco":
        has_input = (getattr(args, "input_dir", None)
                     or getattr(args, "sample_table", None))
        has_fasta_dir = bool(params.get("busco_fasta_dir"))

        if has_fasta_dir:
            fasta_dir = Path(params["busco_fasta_dir"])
            if not fasta_dir.is_dir():
                logging.error(
                    f"[BUSCO] --busco_fasta_dir not found: {fasta_dir}\n"
                    f"  Example: --busco_fasta_dir /my/assemblies/")
                return 1
            fa_list = (sorted(fasta_dir.glob("*.fasta")) +
                       sorted(fasta_dir.glob("*.fa")) +
                       sorted(fasta_dir.glob("*.fna")))
            for fa in fa_list:
                if fa.is_file() and fa.stat().st_size > 0:
                    sid = safe_name(re.sub(r"\.(fasta|fa|fna)$", "", fa.name, flags=re.I))
                    samples.append(Sample(sid=sid, r1="", r2=""))
            if not samples:
                logging.error(f"[BUSCO] no fasta found in {fasta_dir}")
                return 1
            logging.info(f"[BUSCO] fasta_dir={fasta_dir}  found {len(samples)}")

        elif has_input:
            try:
                samples, warns = collect_samples(args, suffix_pairs, require_fastq=False)
            except RuntimeError as e:
                logging.error(str(e))
                return 1
            for w in warns:
                logging.warning(f"[INPUT] {w}")
            valid: List[Sample] = []
            for s in samples:
                if s.r1:
                    valid.append(s)
                else:
                    spades_sub = spades_sample_dir(out_root, params, s.sid)
                    if spades_success_files(spades_sub, params["data_type"]):
                        valid.append(s)
                    else:
                        logging.warning(f"[BUSCO_SKIP] {s.sid}: no r1/spades output")
            samples = valid

        else:
            if params["direct_busco_output"]:
                scan_root = out_root
            else:
                scan_root = out_root / params["spades_dir"]
            if not scan_root.is_dir():
                logging.error(
                    f"[BUSCO] scan dir not found: {scan_root}\n"
                    f"  Provide --input_dir, --sample_table, or --busco_fasta_dir\n"
                    f"  Example: {sys.argv[0]} busco -o {out_root} "
                    f"--busco_fasta_dir /assemblies --busco_lineage bacteria_odb10")
                return 1
            suffix = params.get("suffix_spades_out", "") or ""
            skip_dirs = {out_root / params["fastp_dir"], out_root / params["busco_dir"]}
            # Also build a set of dir names to skip (pipeline dirs + hidden dirs)
            skip_names: Set[str] = {
                params["fastp_dir"], params["busco_dir"],
                ".snakemake", ".snakemake_generated", "__pycache__",
            }
            for d in sorted(scan_root.iterdir()):
                if not d.is_dir() or d in skip_dirs:
                    continue
                dname = d.name
                if dname.startswith(".") or dname in skip_names:
                    continue
                sid = dname[: -len(suffix)] if suffix and dname.endswith(suffix) else dname
                samples.append(Sample(sid=sid, r1="", r2=""))
            logging.info(
                f"[BUSCO_STANDALONE] auto-found {len(samples)} sample dirs "
                f"under {scan_root}")

        if not samples:
            logging.error(
                "[BUSCO] no samples to process.\n"
                "  Hint: use --busco_fasta_dir, -i, or --sample_table")
            return 1

    else:
        try:
            samples, missing = collect_samples(args, suffix_pairs, require_fastq=True)
        except RuntimeError as e:
            logging.error(str(e))
            return 1
        if missing:
            logging.warning(
                f"[INPUT] {len(missing)} samples with missing reads: "
                f"{', '.join(missing[:20])}" + ("..." if len(missing) > 20 else ""))
        if not samples:
            input_dir = getattr(args, "input_dir", None)
            if input_dir:
                input_p = Path(input_dir)
                if input_p.is_dir():
                    fq_files = (list(input_p.glob("*.fastq.gz")) +
                                list(input_p.glob("*.fq.gz")) +
                                list(input_p.glob("**/*.fastq.gz"))[:20])
                    if fq_files:
                        logging.error(
                            f"[INPUT] No valid paired reads detected.\n"
                            f"  Found files like: {fq_files[0].name}\n"
                            f"  but none matched the suffix pattern(s).\n"
                            f"  Try: --suffix_fq auto\n"
                            f"  Or specify: --suffix_fq _R1_001.fastq.gz _R2_001.fastq.gz\n"
                            f"  Presets: {', '.join(sorted(_SUFFIX_PRESETS.keys()))}")
                    else:
                        logging.error(
                            f"[INPUT] No fastq files in {input_dir}\n"
                            f"  If reads are in subdirs: --fq_position subdir")
                else:
                    logging.error(f"[INPUT] --input_dir does not exist: {input_dir}")
            else:
                logging.error(
                    "[INPUT] No input specified.\n"
                    "  Use: -i /path/to/reads  or  --sample_table samples.tsv")
            return 1

    logging.info(f"[PIPELINE] subcommand={args.cmd}  total_samples={len(samples)}")

    # ============================================================
    # v6.3: UNIFIED PRE-FLIGHT (ALL checks before anything runs)
    # ============================================================
    use_smk = getattr(args, "use_snakemake", False)
    auto_yes = getattr(args, "auto_yes", False)
    dryrun = bool(getattr(args, "dryrun", False))

    if not getattr(args, "no_preflight", False):
        rc = preflight_full_check(
            cmd=args.cmd, params=params, samples=samples,
            out_root=out_root, use_smk=use_smk)
        if rc != 0:
            return rc
    else:
        logging.warning("[PREFLIGHT] Skipped (--no_preflight). Not recommended.")

    # ── v6.5: print status table for all samples before execution plan ──
    _status_rows_main = scan_all_status(samples, out_root, params, args.cmd)
    print_status_table(_status_rows_main, args.cmd, params, out_root,
                       max_rows=50)

    # Print execution plan
    print_execution_plan(
        cmd=args.cmd, params=params, samples=samples,
        out_root=out_root, use_smk=use_smk)

    # Ask for confirmation (unless --yes or --background)
    if not auto_yes and not dryrun:
        if not ask_user_confirmation(auto_yes=False):
            return 0
    elif dryrun:
        print("  --dryrun mode: will simulate without executing.\n")
    else:
        print("  --yes flag detected, proceeding automatically.\n")

    # ============================================================
    # v6.3: Pilot run + auto-background
    # ============================================================
    auto_bg = getattr(args, "auto_bg", False)
    pilot_n = max(1, getattr(args, "pilot_n", 5))

    if auto_bg and not dryrun and len(samples) > pilot_n:
        pilot_samples = samples[:pilot_n]
        logging.info(
            f"[AUTO_BG] Pilot run: {pilot_n} samples in foreground first")

        # Run pilot samples
        if use_smk:
            rc_pilot = run_snakemake_pipeline(
                pilot_samples, out_root, params,
                cmd=args.cmd, dryrun=False)
        else:
            rc_pilot = run_direct_pipeline(
                pilot_samples, out_root, params,
                cmd=args.cmd, dryrun=False)

        if rc_pilot != 0:
            logging.error(
                f"[AUTO_BG] Pilot run FAILED (rc={rc_pilot}) on first "
                f"{pilot_n} samples. NOT launching background.\n"
                f"  Fix the issues above and re-run.")
            return rc_pilot

        logging.info(
            f"[AUTO_BG] Pilot run passed ({pilot_n}/{pilot_n}). "
            f"Launching remaining {len(samples) - pilot_n} samples in background...")

        # Re-launch the full pipeline in background (skip-logic will
        # handle the already-done pilot samples automatically).
        return launch_background(out_root)

    # ============================================================
    # v6: dispatch to snakemake or threadpool backend
    # ============================================================
    if use_smk:
        logging.info("[BACKEND] Using Snakemake backend")
        return run_snakemake_pipeline(
            samples, out_root, params,
            cmd=args.cmd, dryrun=dryrun)
    else:
        logging.info("[BACKEND] Using built-in ThreadPoolExecutor backend")
        return run_direct_pipeline(
            samples, out_root, params,
            cmd=args.cmd, dryrun=dryrun)


if __name__ == "__main__":
    raise SystemExit(main())
