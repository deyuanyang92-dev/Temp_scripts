#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
blast2_metadatav2.12.py

Batch BLAST (local/online) -> optional TopN -> optional Entrez metadata -> merge summary.

CODEX PATCH NOTES (2026-04-22)
==============================
Codex 修改: 在不分叉新脚本的前提下继续维护 v2.12，使其覆盖 v2.1/v2.7/v2.8
仍有价值的功能，并在修改处保留“原代码/原逻辑”注释摘要。
理由:
  1) v2.7 的 metadata-only accession-list 读取能力在 v2.12 中缺失；
  2) v2.12 的 topN-only、blast2topN、失败退出码、dry-run fallback、CSV/TSV
     后备输出仍有运行级 bug；
  3) v2.12 已包含旧版的本地 metadata、伪 accession、Snakemake、限速和
     outfmt 修复，是最适合作为唯一主脚本的基线。

CHANGES FROM v2.11 (code-review round 3 - final production fixes)
==================================================================
Fix-G  [Critical] BLAST vs. Entrez rate-limit collision:  Split NCBI_LIMITER
         into ENTREZ_LIMITER (respects API key: 3 or 10 req/sec) and 
         BLAST_LIMITER (rigid 1 req / 3 sec, no API key support).  The NCBI
         BLAST API ([https://blast.ncbi.nlm.nih.gov/Blast.cgi](https://blast.ncbi.nlm.nih.gov/Blast.cgi)) is a separate 
         backend that does NOT honor Entrez API keys.  Sending BLAST requests
         at 10/sec with an API key triggers immediate HTTP 429/503 errors and
         permanent IP bans.

Fix-H  [Critical] csv_output always written (4 broken code paths fixed):
         • blast-only mode: now reads the merged BLAST file and writes it to 
           csv_output before returning — previously wrote nothing, causing 
           Snakemake's `test -s {output.csv}` to fail and --resume to never 
           trigger.
         • topN-only / metadata-only / blast2*/no-pandas: all wrote to 
           csv_output.with_suffix(".tsv") — now write directly to csv_output,
           fixing the same Snakemake sentinel failure and resume-detection bug.

Fix-I  [Minor] Stale double-comment block above _RateLimiter cleaned up
         (leftover fragment from v2.11 edit merged into clean single block).

Fix-J  [Minor] --help epilog examples still referenced 
         blast2_metadatav2.10_fixed.py — updated to blast2_metadatav2.12.py.

CHANGES FROM v2.10 (code-review round 2)
=========================================
Fix-A  [Critical] Inter-process rate limiter:  _RateLimiter is now re-initialised
         inside main() AFTER args are parsed.  When running under Snakemake,
         `_smk_active_jobs` is injected into child_args.json so every child
         process divides the global NCBI budget (3 or 10 req/sec) evenly across
         all concurrently running jobs — preventing IP bans.

Fix-B  [Critical] Snakemake CLI clash:  The redundant `-j` flag has been
         removed from the snakemake command.  Modern Snakemake (v8+) treats
         `--cores` / `-c` and `--jobs` / `-j` as strict aliases; passing both
         causes a crash.  Only `--cores` is now passed.

Fix-C  [Critical] Biopython XML parser resource leak:  `NCBIXML.parse(handle)`
         is now wrapped in `try...finally: handle.close()` so the socket is
         always closed even when NCBI returns an HTML error page that triggers
         an ExpatError during XML parsing.

Fix-D  [Minor] API-key rate-limit under-utilisation:  When a valid `--api_key`
         is supplied, NCBI allows 10 req/sec instead of 3.  The limiter now
         reflects this automatically.

Fix-E  [Minor] Entrez efetch handle leak in `do_fetch`:  `handle.close()` is
         now inside a `try...finally` block so it is always called even when
         `Entrez.read()` raises.

Fix-F  [Minor] Removed unused `field` import from `dataclasses`.

CHANGES FROM v2.9 (code-review round 1)
=========================================
Fix-1  [Critical-A] Snakemake child-process arg passthrough via --args_json.
Fix-2  [Critical-B] _parse_outfmt6_line dynamic column list.
Fix-3  [Critical-C] build_outfmt_string returns (fmt_string, col_names).
Fix-4  [Critical-D/E] Global thread-safe _RateLimiter for all NCBI calls.
Fix-5  [Design-A] Snakefile template updated to use --args_json.
Fix-6  [Design-D] col_names propagated through parse/write pipeline.
Fix-7  Minor build_outfmt_string validation kept; helpers untouched.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import glob
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import time
import threading
import shlex
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

# Optional pandas
try:
    import pandas as pd  # type: ignore
    HAVE_PANDAS = True
except Exception:
    pd = None  # type: ignore
    HAVE_PANDAS = False

# Optional psutil
try:
    import psutil  # type: ignore
    HAVE_PSUTIL = True
except Exception:
    psutil = None  # type: ignore
    HAVE_PSUTIL = False

# Biopython required
try:
    from Bio import Entrez, SeqIO  # type: ignore
    from Bio.Blast import NCBIWWW, NCBIXML  # type: ignore
except Exception as e:
    raise SystemExit(
        "ERROR: Biopython is required. Install with: pip install biopython\n"
        f"Original error: {e}"
    )

# ─────────────────────────────────────────────
# logging
# ─────────────────────────────────────────────
logger = logging.getLogger("blast2_metadata")
logger.setLevel(logging.INFO)
_handler = logging.StreamHandler(sys.stdout)
_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logger.handlers.clear()
logger.addHandler(_handler)

BLAST_HEADERS: List[str] = [
    "qseqid", "sseqid", "stitle", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs",
]

NT_META_COLS: List[str] = [
    "accession", "version", "locus", "length", "definition",
    "molecule_type", "topology", "division", "update_date", "create_date",
    "organism", "taxonomy",
    "reference_titles", "reference_authors", "reference_journals", "comment",
]

NUMERIC_FLOAT_COLS = {"pident", "evalue", "qcovs", "bitscore"}
NUMERIC_INT_COLS   = {"length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send"}

# Codex 修改: 增加集中式扩展名、模式和 accession DB 识别常量。
# 理由: 修复 metadata-only/topN-only 输入发现不完整，以及 Entrez 蛋白/核酸库误判问题。
FASTA_EXTENSIONS_FALLBACK = [".fasta", ".fas", ".faa", ".fa", ".fna"]
ACCESSION_LIST_EXTENSIONS = [".txt", ".tsv", ".csv", ".list", ".acc", ".accessions"]
BLAST_TABLE_EXTENSIONS    = [".tsv", ".tab", ".txt", ".csv"]
BLAST_REQUIRED_MODES      = {"blast-only", "blast2topN", "blast2meta", "blast2topN2metadata"}

PROTEIN_ACCESSION_PREFIXES = (
    "NP_", "YP_", "WP_", "XP_", "AP_", "ZP_", "CAA", "AAB", "AAC", "AAD",
    "AAH", "AAM", "AAN", "AAO", "AAP", "AAQ", "AAR", "AAS", "AAT",
)
NUCLEOTIDE_ACCESSION_PREFIXES = (
    "NC_", "NG_", "NM_", "NR_", "NT_", "NW_", "NZ_", "XM_", "XR_", "AC_",
    "AP", "CP", "CM", "JX", "KC", "KF", "KP", "KU", "KY", "MK", "MN",
    "MT", "MW", "MZ", "OP", "OQ", "OR", "PQ",
)


# ─────────────────────────────────────────────
# Fix-G: Split NCBI rate limiter into ENTREZ and BLAST
# ─────────────────────────────────────────────
# The NCBI API Key system *only applies to Entrez E-utilities* (esearch, efetch).
# The NCBI BLAST API (https://blast.ncbi.nlm.nih.gov/Blast.cgi) is a completely
# separate backend that does NOT support API keys and enforces a strict limit
# of 1 request every 3 seconds. Sharing a single rate limiter between them
# causes immediate IP bans when users supply an API key (bumping the limit to
# 10/sec) and fire BLAST requests at that rate.
#
# Both limiters are RE-INITIALISED inside main() once real args are known:
#   • api_key → bumps ENTREZ limit to 10 req/sec (Fix-D)
#   • Snakemake child processes divide budgets by _smk_active_jobs (Fix-A)
class _RateLimiter:
    """Thread-safe token-bucket rate limiter."""

    def __init__(self, max_calls: int, period: float) -> None:
        self._max = max_calls
        self._period = period
        self._lock = threading.Lock()
        self._timestamps: List[float] = []

    def acquire(self) -> None:
        while True:
            with self._lock:
                now = time.monotonic()
                # purge calls older than one period
                self._timestamps = [t for t in self._timestamps if now - t < self._period]
                if len(self._timestamps) < self._max:
                    self._timestamps.append(now)
                    return
                oldest = self._timestamps[0]
                wait = self._period - (now - oldest)
            time.sleep(max(0.05, wait))


def _make_entrez_limiter(api_key: str = "", active_jobs: int = 1) -> "_RateLimiter":
    """
    Build Entrez E-utilities rate limiter.

    Parameters
    ----------
    api_key     : non-empty string → NCBI allows 10 req/sec instead of 3
    active_jobs : number of concurrently running Snakemake child processes;
                  the per-process share of the global budget is divided evenly.
    """
    max_reqs    = 10 if api_key else 3
    active_jobs = max(1, int(active_jobs))
    # Give each process a fair share of the total NCBI budget, with a tiny
    # extra margin (1.05×) to absorb clock jitter.
    safe_period = 1.05 * active_jobs
    return _RateLimiter(max_calls=max_reqs, period=safe_period)


def _make_blast_limiter(active_jobs: int = 1) -> "_RateLimiter":
    """
    Build NCBI BLAST API rate limiter.
    
    CRITICAL: The BLAST API enforces a STRICT limit of 1 request every 3 seconds,
    regardless of API keys. To be safe across distributed clock jitter, we use 
    3.1 seconds per request.
    
    Parameters
    ----------
    active_jobs : number of concurrently running Snakemake child processes;
                  each gets 1 request per (3.1 × active_jobs) seconds.
    """
    active_jobs = max(1, int(active_jobs))
    # STRICT limit: 1 request every 3 seconds. 
    # To be safe across distributed clock jitter, we use 3.1 seconds.
    safe_period = 3.1 * active_jobs
    return _RateLimiter(max_calls=1, period=safe_period)


# Module-level sentinels — replaced in main() after args are parsed.
ENTREZ_LIMITER: _RateLimiter = _make_entrez_limiter()
BLAST_LIMITER: _RateLimiter = _make_blast_limiter()


# ─────────────────────────────────────────────
# generic helpers
# ─────────────────────────────────────────────
def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def which(exe: str) -> Optional[str]:
    try:
        return shutil.which(exe)
    except Exception:
        return None


def is_file_complete(file_path: Union[str, Path]) -> bool:
    p = Path(file_path)
    return p.exists() and p.is_file() and p.stat().st_size > 0


def _shlex_split(s: str) -> List[str]:
    try:
        return shlex.split(s)
    except Exception:
        return s.split()


def _extra_has_outfmt(extra: str) -> bool:
    if not extra:
        return False
    toks = _shlex_split(extra)
    for t in toks:
        if t.strip() in {"-outfmt", "--outfmt"}:
            return True
    return any(t.startswith("-outfmt") for t in toks)


def run_cmd(
    cmd: List[str],
    cwd: Optional[Path] = None,
    retries: int = 1,
    retry_sleep: int = 10,
    env: Optional[Dict[str, str]] = None,
) -> Tuple[int, str, str]:
    last = (1, "", "")
    sleep = max(1, int(retry_sleep))
    for attempt in range(1, retries + 1):
        logger.info(f"[CMD] {' '.join(cmd)}")
        p = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            errors="replace",
            env=env,
        )
        last = (p.returncode, p.stdout or "", p.stderr or "")
        if p.returncode == 0:
            return last
        logger.warning(f"[CMD_FAIL] attempt {attempt}/{retries} rc={p.returncode}")
        if attempt < retries:
            time.sleep(sleep)
            sleep = int(sleep * 1.8)
    return last


# ─────────────────────────────────────────────
# Fix-3: build_outfmt_string → returns (str, List[str])
# ─────────────────────────────────────────────
def build_outfmt_string(outfmt_value: str) -> Tuple[str, List[str]]:
    """
    Returns (blast_outfmt_arg, list_of_column_names).

    - "6" / "7" / ""  → use all BLAST_HEADERS
    - "6 col1 col2 …" → validate required headers present, use provided list
    - anything else   → pass through unchanged, assume BLAST_HEADERS order
    """
    s = (outfmt_value or "").strip()

    # Fix-3: default empty / bare digit to standard 14-column format
    if s in {"", "6", "7"}:
        fmt_num = s if s in {"6", "7"} else "6"
        fmt_str = fmt_num + " " + " ".join(BLAST_HEADERS)
        return fmt_str, list(BLAST_HEADERS)

    toks = s.split()
    if toks and toks[0] in {"6", "7"}:
        fields = toks[1:]
        if not fields:
            # "6" or "7" with no columns → use defaults
            fmt_str = toks[0] + " " + " ".join(BLAST_HEADERS)
            return fmt_str, list(BLAST_HEADERS)
        missing = [h for h in BLAST_HEADERS if h not in fields]
        if missing:
            raise SystemExit(
                f"ERROR: --outfmt_blast '{outfmt_value}' is missing required fields: "
                f"{missing}. Use --outfmt_blast 6 to get all required columns automatically."
            )
        return s, list(fields)

    # Codex 修改: 非 6/7 的 BLAST 输出格式直接报错。
    # 理由: 这个 pipeline 后续解析依赖 tabular 输出；原代码会放行 XML/pairwise，
    #       然后在 TopN/metadata 阶段静默解析为空或错列。
    # 原代码:
    #   # Non-tabular format (0, 5, etc.) – pass through, assume BLAST_HEADERS for parsing
    #   return s, list(BLAST_HEADERS)
    raise SystemExit(
        f"ERROR: --outfmt_blast '{outfmt_value}' is not supported by this pipeline. "
        "Use tabular BLAST output format 6 or 7."
    )


def detect_sequence_type(fasta_file: Union[str, Path], max_records: int = 2) -> str:
    fasta_file = str(fasta_file)
    try:
        with open(fasta_file, "r", encoding="utf-8", errors="ignore") as f:
            it = SeqIO.parse(f, "fasta")
            checked = 0
            for record in it:
                seq = str(record.seq).upper()
                if not seq:
                    continue
                checked += 1
                dna_count = sum(seq.count(x) for x in ["A", "T", "G", "C", "U", "N"])
                if dna_count / max(1, len(seq)) < 0.85:
                    return "protein"
                if checked >= max_records:
                    break
    except Exception:
        pass
    return "nucleotide"


def get_blast_program(seq_type: str) -> str:
    return "blastn" if seq_type == "nucleotide" else "blastp"


def get_blast_db(seq_type: str, args: argparse.Namespace) -> str:
    if args.blast_db != "auto":
        return args.blast_db
    if args.blast_method == "local":
        return args.blast_db_nt if seq_type == "nucleotide" else args.blast_db_nr
    return "nt" if seq_type == "nucleotide" else "nr"


def _to_float(x: object) -> Optional[float]:
    try:
        return float(x)  # type: ignore
    except Exception:
        return None


def _to_int(x: object) -> Optional[int]:
    try:
        return int(float(x))  # type: ignore
    except Exception:
        return None


def safe_int(x: Optional[Union[str, int]], default: int) -> int:
    try:
        return int(x) if x is not None else int(default)
    except Exception:
        return int(default)


def resolve_total_threads(args: argparse.Namespace) -> int:
    cpu_count = max(1, os.cpu_count() or 1)
    requested = safe_int(getattr(args, "threads", 0), 0)
    if requested <= 0:
        return cpu_count
    return max(1, min(requested, cpu_count))


@dataclass
class HostState:
    cpu_count: int
    mem_total_mb: int
    mem_available_mb: int
    cpu_percent: float
    loadavg1: float
    loadavg5: float
    loadavg15: float
    load_ratio_1m: float
    mem_pressure: float
    busy_score: float


def detect_host_state() -> HostState:
    cpu_count = max(1, os.cpu_count() or 1)
    mem_total_mb = 0
    mem_available_mb = 0
    cpu_percent = 0.0

    if HAVE_PSUTIL:
        try:
            vm = psutil.virtual_memory()
            mem_total_mb     = max(1, int(vm.total     / 1024 / 1024))
            mem_available_mb = max(1, int(vm.available / 1024 / 1024))
            cpu_percent      = float(psutil.cpu_percent(interval=0.2))
        except Exception:
            pass

    if mem_total_mb <= 0:
        try:
            if Path("/proc/meminfo").exists():
                with open("/proc/meminfo", "r", encoding="utf-8", errors="ignore") as fh:
                    for line in fh:
                        if line.startswith("MemTotal:"):
                            mem_total_mb = max(1, int(line.split()[1]) // 1024)
                        elif line.startswith("MemAvailable:"):
                            mem_available_mb = max(1, int(line.split()[1]) // 1024)
        except Exception:
            pass

    if mem_total_mb <= 0:
        try:
            pages     = os.sysconf("SC_PHYS_PAGES")
            page_size = os.sysconf("SC_PAGE_SIZE")
            mem_total_mb = max(1, int((pages * page_size) / 1024 / 1024))
        except Exception:
            mem_total_mb = 32768

    if mem_available_mb <= 0:
        mem_available_mb = max(1, int(mem_total_mb * 0.5))

    try:
        la1, la5, la15 = os.getloadavg()
    except Exception:
        la1, la5, la15 = 0.0, 0.0, 0.0

    load_ratio_1m = float(la1) / max(1, cpu_count)
    mem_pressure  = 1.0 - (float(mem_available_mb) / max(1.0, float(mem_total_mb)))

    cpu_ratio  = max(0.0, min(1.5, cpu_percent / 100.0))
    load_ratio = max(0.0, min(1.5, load_ratio_1m))
    mem_ratio  = max(0.0, min(1.5, mem_pressure))
    busy_score = max(cpu_ratio, load_ratio, mem_ratio)

    return HostState(
        cpu_count=cpu_count,
        mem_total_mb=mem_total_mb,
        mem_available_mb=mem_available_mb,
        cpu_percent=cpu_percent,
        loadavg1=float(la1),
        loadavg5=float(la5),
        loadavg15=float(la15),
        load_ratio_1m=load_ratio_1m,
        mem_pressure=mem_pressure,
        busy_score=busy_score,
    )

# ─────────────────────────────────────────────
# accession parsing
# ─────────────────────────────────────────────
def parse_accession_from_sseqid(sseqid: str) -> str:
    s = (sseqid or "").strip()
    if not s:
        return ""
    patterns = [
        r"^ref\|([^|]+)\|",
        r"\|ref\|([^|]+)\|",
        r"^gb\|([^|]+)\|",
        r"\bgb\|([^|]+)\|",
        r"^emb\|([^|]+)\|",
        r"^dbj\|([^|]+)\|",
        r"^sp\|([^|]+)\|",
        r"^tr\|([^|]+)\|",
    ]
    for pat in patterns:
        m = re.search(pat, s)
        if m:
            return m.group(1).strip()

    toks = s.split("|")
    if len(toks) >= 2 and toks[0] in {"ref", "gb", "emb", "dbj", "sp", "tr"}:
        return toks[1].strip()
    if len(toks) >= 4 and toks[0] == "gi" and toks[2] in {"ref", "gb", "emb", "dbj"}:
        return toks[3].strip()
    return s.split()[0].strip()


def accession_base(acc: str) -> str:
    return (acc or "").strip().split(".")[0]


# Codex 修改: 新增输入类型、accession-list 读取和 Entrez DB 猜测工具函数。
# 理由: 旧版 v2.7 有 metadata-only 的 accession 文件读取能力；v2.12 遗漏后导致该模式空跑。
def path_has_extension(path: Union[str, Path], extensions: Sequence[str]) -> bool:
    name = str(path).lower()
    return any(name.endswith(ext.lower()) for ext in extensions)


def looks_like_fasta_path(path: Union[str, Path], extensions: Optional[Sequence[str]] = None) -> bool:
    return path_has_extension(path, extensions or FASTA_EXTENSIONS_FALLBACK)


def strip_blast_suffix(stem: str) -> str:
    for suffix in ["_blast", ".blast", "-blast"]:
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def read_accession_list(file_path: Union[str, Path]) -> List[str]:
    """Read accession IDs from plain text, TSV or CSV accession-list files."""
    p = Path(file_path)
    if not p.exists() or p.stat().st_size == 0:
        return []

    accessions: List[str] = []
    seen = set()

    with open(p, "r", encoding="utf-8", errors="ignore", newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        except Exception:
            dialect = csv.excel_tab if "\t" in sample else csv.excel

        has_delimiter = ("\t" in sample) or ("," in sample)
        if has_delimiter:
            reader = csv.reader(f, dialect)
            header_checked = False
            acc_col: Optional[int] = None
            for row in reader:
                vals = [x.strip() for x in row]
                if not any(vals):
                    continue
                if not header_checked:
                    lowered = [x.lower() for x in vals]
                    for candidate in ["accession", "accession_number", "acc", "id", "sseqid"]:
                        if candidate in lowered:
                            acc_col = lowered.index(candidate)
                            header_checked = True
                            break
                    if acc_col is None:
                        header_checked = True
                    if acc_col is not None and any(x in lowered for x in ["accession", "sseqid", "id"]):
                        continue
                raw = vals[acc_col] if acc_col is not None and acc_col < len(vals) else vals[0]
                acc = accession_base(parse_accession_from_sseqid(raw))
                if acc and acc not in seen:
                    seen.add(acc)
                    accessions.append(acc)
        else:
            for line in f:
                raw = line.strip()
                if not raw or raw.startswith("#"):
                    continue
                acc = accession_base(parse_accession_from_sseqid(raw.split()[0]))
                if acc and acc not in seen:
                    seen.add(acc)
                    accessions.append(acc)

    return accessions


def guess_entrez_db_for_accession(acc: str, default_db: str = "nucleotide") -> str:
    base = accession_base(acc).upper()
    if base.startswith(PROTEIN_ACCESSION_PREFIXES):
        return "protein"
    if base.startswith(NUCLEOTIDE_ACCESSION_PREFIXES):
        return "nucleotide"
    return default_db


def is_pseudo_accession(acc_base: str, prefixes: List[str]) -> bool:
    if not acc_base:
        return False
    for p in prefixes:
        p = p.strip()
        if p and acc_base.startswith(p):
            return True
    return False


def is_local_or_pseudo_accession(
    acc: str, local_meta: Dict[str, Dict[str, str]], prefixes: List[str]
) -> bool:
    if not acc:
        return False
    acc_full = (acc or "").strip()
    acc_base = accession_base(acc_full)
    if acc_full in local_meta or acc_base in local_meta:
        return True
    return is_pseudo_accession(acc_base, prefixes)


# ─────────────────────────────────────────────
# Fix-2: BLAST output parsing with dynamic col list
# ─────────────────────────────────────────────
def _parse_outfmt6_line(
    parts: List[str],
    col_names: List[str],
) -> Optional[Dict[str, object]]:
    """
    Parse one tab-split line according to *col_names* (the actual column list
    returned by build_outfmt_string).

    Handles three cases robustly:
      • n == expected  → direct 1-to-1 mapping
      • n > expected and 'stitle' in col_names → extra tabs inside stitle; merge middle
      • n > expected and no stitle             → extra custom columns appended; take first
    """
    n        = len(parts)
    expected = len(col_names)

    if n < expected:
        return None  # malformed line

    if n == expected:
        cols = parts
    elif "stitle" in col_names:
        stitle_idx = col_names.index("stitle")
        post_count = expected - stitle_idx - 1
        pre  = parts[:stitle_idx]
        post = parts[n - post_count:] if post_count > 0 else []
        stitle_merged = "\t".join(parts[stitle_idx: n - post_count if post_count else n])
        cols = pre + [stitle_merged] + post
    else:
        # extra columns at end (user added qlen, slen, etc.); keep first `expected`
        cols = parts[:expected]

    r: Dict[str, object] = {}
    for i, k in enumerate(col_names):
        v = cols[i] if i < len(cols) else None
        if k in NUMERIC_FLOAT_COLS:
            r[k] = _to_float(v)
        elif k in NUMERIC_INT_COLS:
            r[k] = _to_int(v)
        else:
            r[k] = v
    return r


def parse_blast_results_py(
    blast_output: Union[str, Path],
    col_names: Optional[List[str]] = None,
) -> List[Dict[str, object]]:
    """Parse a tabular BLAST output file.  col_names defaults to BLAST_HEADERS."""
    if col_names is None:
        col_names = BLAST_HEADERS
    p = Path(blast_output)
    if not p.exists() or p.stat().st_size == 0:
        return []
    rows: List[Dict[str, object]] = []
    with open(p, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            r = _parse_outfmt6_line(parts, col_names)
            if r is None:
                continue
            rows.append(r)
    return rows


# Codex 修改: default_headers 可由调用方指定。
# 理由: 原函数空结果时总是 BLAST_HEADERS，metadata-only 空结果会写出 BLAST 表头。
def rows_to_df(
    rows: List[Dict[str, object]], default_headers: Optional[Sequence[str]] = None
) -> "pd.DataFrame":
    assert HAVE_PANDAS and pd is not None
    if not rows:
        return pd.DataFrame(columns=list(default_headers or BLAST_HEADERS))
    return pd.DataFrame(rows)


def merge_blast_tsv_text(blast_files: Sequence[Union[str, Path]], merged_out: Path) -> None:
    ensure_dir(merged_out.parent)
    with open(merged_out, "w", encoding="utf-8") as out:
        for fp in blast_files:
            p = Path(fp)
            if not p.exists() or p.stat().st_size == 0:
                continue
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if not line.strip() or line.startswith("#"):
                        continue
                    out.write(line)


# ─────────────────────────────────────────────
# TopN
# ─────────────────────────────────────────────
def get_topN_hits_py(
    rows: List[Dict[str, object]], topN: int, sort_by: str
) -> List[Dict[str, object]]:
    if not rows:
        return rows

    def fnum(v, default: float) -> float:
        if v is None:
            return default
        try:
            return float(v)
        except Exception:
            return default

    def key_identity(r):
        return (str(r.get("qseqid", "")), -fnum(r.get("pident"), -1e9),
                -fnum(r.get("qcovs"), -1e9), fnum(r.get("evalue"), 1e9))

    def key_evalue(r):
        return (str(r.get("qseqid", "")), fnum(r.get("evalue"), 1e9),
                -fnum(r.get("pident"), -1e9), -fnum(r.get("qcovs"), -1e9))

    def key_cov(r):
        return (str(r.get("qseqid", "")), -fnum(r.get("qcovs"), -1e9),
                -fnum(r.get("pident"), -1e9), fnum(r.get("evalue"), 1e9))

    if sort_by == "evalue":
        sorted_rows = sorted(rows, key=key_evalue)
    elif sort_by == "coverage":
        sorted_rows = sorted(rows, key=key_cov)
    else:
        sorted_rows = sorted(rows, key=key_identity)

    out: List[Dict[str, object]] = []
    seen: Dict[str, int] = {}
    for r in sorted_rows:
        q = str(r.get("qseqid", ""))
        c = seen.get(q, 0)
        if c < topN:
            out.append(r)
            seen[q] = c + 1
    return out


def extract_accession_numbers_from_rows(rows: List[Dict[str, object]]) -> List[str]:
    clean: List[str] = []
    for r in rows:
        acc  = parse_accession_from_sseqid(str(r.get("sseqid", "")))
        a    = accession_base(acc)
        if a and a not in clean:
            clean.append(a)
    return clean


# ─────────────────────────────────────────────
# FASTA batching
# ─────────────────────────────────────────────
@dataclass
class BatchInfo:
    batch_id: str
    fasta_path: Path


def iter_fasta_batches(
    fasta_file: Path, temp_dir: Path, batch_size: int
) -> List[BatchInfo]:
    ensure_dir(temp_dir)
    batches: List[BatchInfo] = []
    bidx = 0
    cur_records = []

    def flush():
        nonlocal bidx, cur_records
        if not cur_records:
            return
        bid     = f"batch_{bidx:05d}"
        out_fa  = temp_dir / f"{bid}__{fasta_file.name}"
        with open(out_fa, "w", encoding="utf-8") as fh:
            SeqIO.write(cur_records, fh, "fasta")
        batches.append(BatchInfo(batch_id=bid, fasta_path=out_fa))
        bidx += 1
        cur_records = []

    with open(fasta_file, "r", encoding="utf-8", errors="ignore") as f:
        for rec in SeqIO.parse(f, "fasta"):
            cur_records.append(rec)
            if len(cur_records) >= batch_size:
                flush()
        flush()

    return batches


# ─────────────────────────────────────────────
# BLAST runners
# ─────────────────────────────────────────────
def compute_parallel_params(args: argparse.Namespace, mode: str) -> Tuple[int, int]:
    total_threads = max(1, int(args.threads))
    pj = int(args.parallel_jobs) if args.parallel_jobs is not None else 0

    if mode == "online":
        if pj <= 0:
            pj = max(1, int(args.online_max_concurrent))
        pj = max(1, min(pj, int(args.online_max_concurrent)))
        return pj, 1

    if pj <= 0:
        pj = min(4, total_threads)
    pj = max(1, min(pj, total_threads))
    if args.local_threads_per_job and int(args.local_threads_per_job) > 0:
        tpj = int(args.local_threads_per_job)
    else:
        tpj = int(math.ceil(total_threads / pj))
    tpj = max(1, tpj)
    return pj, tpj


def run_local_blast_one_batch(
    batch: BatchInfo,
    seq_type: str,
    args: argparse.Namespace,
    blast_dir: Path,
    threads_per_job: int,
    col_names: List[str],
) -> Tuple[str, bool, Path, List[str]]:
    blast_dir = ensure_dir(blast_dir)
    out_file  = blast_dir / f"{batch.batch_id}_blast.tsv"
    if args.resume and is_file_complete(out_file):
        return (batch.batch_id, True, out_file, col_names)

    program = args.blast_exe if args.blast_exe != "auto" else get_blast_program(seq_type)
    if not which(program):
        raise FileNotFoundError(f"BLAST program not found in PATH: {program}")

    db = get_blast_db(seq_type, args)
    # Fix-3: unpack tuple
    outfmt_str, col_names = build_outfmt_string(args.outfmt_blast)

    if _extra_has_outfmt(args.extra_para2_blast):
        raise SystemExit(
            "ERROR: --extra_para2_blast must NOT include -outfmt. "
            "Remove it and use --outfmt_blast instead."
        )

    cmd = [
        program,
        "-query", str(batch.fasta_path),
        "-db",    str(db),
        "-outfmt", outfmt_str,
        "-out",   str(out_file),
        "-evalue", str(args.evalue),
        "-max_target_seqs", str(args.max_target_seqs),
        "-num_threads", str(max(1, int(threads_per_job))),
    ]
    if args.extra_para2_blast:
        cmd += _shlex_split(args.extra_para2_blast)

    rc, so, se = run_cmd(cmd, retries=args.retries, retry_sleep=args.retry_sleep)
    if rc != 0:
        logger.error(f"[LOCAL_BLAST_FAIL] {batch.batch_id}\n{se[-2000:]}")
        return (batch.batch_id, False, out_file, col_names)
    return (batch.batch_id, True, out_file, col_names)


def run_remote_blast_one_batch(
    batch: BatchInfo,
    seq_type: str,
    args: argparse.Namespace,
    blast_dir: Path,
    col_names: List[str],
) -> Tuple[str, bool, Path, List[str]]:
    blast_dir = ensure_dir(blast_dir)
    out_file  = blast_dir / f"{batch.batch_id}_blast.tsv"
    if args.resume and is_file_complete(out_file):
        return (batch.batch_id, True, out_file, col_names)

    program = args.blast_exe if args.blast_exe != "auto" else get_blast_program(seq_type)
    if not which(program):
        raise FileNotFoundError(
            f"BLAST+ program not found in PATH: {program}\n"
            "Install BLAST+ (blastn/blastp) or use --blast_method local."
        )

    db = get_blast_db(seq_type, args)
    outfmt_str, col_names = build_outfmt_string(args.outfmt_blast)

    if _extra_has_outfmt(args.extra_para2_blast):
        raise SystemExit(
            "ERROR: --extra_para2_blast must NOT include -outfmt. "
            "Remove it and use --outfmt_blast instead."
        )

    cmd = [
        program, "-remote",
        "-query", str(batch.fasta_path),
        "-db",    str(db),
        "-outfmt", outfmt_str,
        "-out",   str(out_file),
        "-evalue", str(args.evalue),
        "-max_target_seqs", str(args.max_target_seqs),
    ]
    if args.extra_para2_blast:
        cmd += _shlex_split(args.extra_para2_blast)

    # Fix-G: use BLAST_LIMITER for remote BLAST calls
    BLAST_LIMITER.acquire()
    rc, so, se = run_cmd(cmd, retries=args.retries, retry_sleep=args.retry_sleep)
    if rc != 0:
        logger.error(f"[REMOTE_BLAST_FAIL] {batch.batch_id}\n{se[-2000:]}")
        return (batch.batch_id, False, out_file, col_names)

    if args.online_delay_sec > 0:
        time.sleep(float(args.online_delay_sec))
    return (batch.batch_id, True, out_file, col_names)


def _parse_extra_params_for_qblast(extra: str) -> Dict[str, Union[str, bool]]:
    params: Dict[str, Union[str, bool]] = {}
    if not extra:
        return params
    xs = _shlex_split(extra)
    i = 0
    while i < len(xs):
        if xs[i].startswith("-"):
            k = xs[i].lstrip("-")
            if i + 1 < len(xs) and not xs[i + 1].startswith("-"):
                params[k] = xs[i + 1]
                i += 2
            else:
                params[k] = True
                i += 1
        else:
            i += 1
    return params


def entrez_setup(args: argparse.Namespace) -> None:
    if args.email:
        Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
    if args.entrez_tool:
        Entrez.tool = args.entrez_tool


def run_biopython_qblast_one_batch(
    batch: BatchInfo,
    seq_type: str,
    args: argparse.Namespace,
    blast_dir: Path,
    col_names: List[str],
) -> Tuple[str, bool, Path, List[str]]:
    blast_dir = ensure_dir(blast_dir)
    out_file  = blast_dir / f"{batch.batch_id}_blast.tsv"
    if args.resume and is_file_complete(out_file):
        return (batch.batch_id, True, out_file, col_names)

    if _extra_has_outfmt(args.extra_para2_blast):
        raise SystemExit(
            "ERROR: --extra_para2_blast must NOT include -outfmt. "
            "Remove it and use --outfmt_blast instead."
        )

    entrez_setup(args)
    qblast_program = get_blast_program(seq_type)
    db    = get_blast_db(seq_type, args)
    extra = _parse_extra_params_for_qblast(args.extra_para2_blast)

    seqs = list(SeqIO.parse(str(batch.fasta_path), "fasta"))
    with open(out_file, "w", encoding="utf-8") as out:
        for idx, rec in enumerate(seqs, start=1):
            qseq = str(rec.seq)
            qid  = rec.id
            qlen = len(rec.seq)

            delay = 5
            ok    = False
            for attempt in range(1, args.retries + 1):
                try:
                    logger.info(
                        f"[QBLAST] {batch.batch_id} {idx}/{len(seqs)} {qid} attempt {attempt}"
                    )
                    # Fix-G: use BLAST_LIMITER for qblast
                    # Fix-C: always close handle even if XML parsing fails
                    BLAST_LIMITER.acquire()
                    handle = NCBIWWW.qblast(
                        program=qblast_program,
                        database=db,
                        sequence=qseq,
                        expect=args.evalue,
                        hitlist_size=args.max_target_seqs,
                        **extra,
                    )
                    try:
                        blast_records = list(NCBIXML.parse(handle))
                    finally:
                        handle.close()

                    for br in blast_records:
                        for aln in br.alignments:
                            for hsp in aln.hsps:
                                qcov   = (hsp.query_end - hsp.query_start + 1) / max(1, qlen) * 100.0
                                subj_id = getattr(aln, "accession", "") or getattr(aln, "hit_id", "")
                                title   = getattr(aln, "title", "") or ""
                                pid     = (hsp.identities / max(1, hsp.align_length)) * 100.0
                                mism    = hsp.align_length - hsp.identities
                                gaps    = getattr(hsp, "gaps", 0) or 0
                                out.write(
                                    f"{qid}\t{subj_id}\t{title}\t{pid:.2f}\t{hsp.align_length}\t"
                                    f"{mism}\t{gaps}\t{hsp.query_start}\t{hsp.query_end}\t"
                                    f"{hsp.sbjct_start}\t{hsp.sbjct_end}\t"
                                    f"{hsp.expect}\t{hsp.score}\t{qcov:.2f}\n"
                                )
                    ok = True
                    break
                except Exception as e:
                    logger.warning(f"[QBLAST_FAIL] {qid} attempt {attempt}: {e}")
                    if attempt < args.retries:
                        time.sleep(delay)
                        delay = int(delay * 1.8)

            if not ok:
                return (batch.batch_id, False, out_file, col_names)

            if args.online_delay_sec > 0:
                time.sleep(float(args.online_delay_sec))

    return (batch.batch_id, True, out_file, col_names)



# ─────────────────────────────────────────────
# local meta loader
# ─────────────────────────────────────────────
def _add_meta_alias(meta: Dict[str, Dict[str, str]], key: str, row: Dict[str, str]) -> None:
    k = (key or "").strip()
    if not k:
        return
    meta[k] = row
    kb = accession_base(k)
    if kb:
        meta.setdefault(kb, row)


def load_local_meta(
    meta_tsv: Union[str, Path], key_col: str
) -> Tuple[Dict[str, Dict[str, str]], str, List[str]]:
    p = Path(meta_tsv)
    if not p.exists() or p.stat().st_size == 0:
        return {}, key_col, []

    # Warn for large files
    size_mb = p.stat().st_size / 1024 / 1024
    if size_mb > 500:
        logger.warning(
            f"[LOCAL_META] {meta_tsv} is {size_mb:.0f} MB. "
            "Each parallel job loads its own copy. Consider splitting or converting to SQLite."
        )

    meta: Dict[str, Dict[str, str]] = {}
    key_col_use = key_col
    fieldnames:  List[str] = []

    with open(p, "r", encoding="utf-8", errors="ignore") as f:
        r = csv.DictReader(f, delimiter="\t")
        fieldnames = list(r.fieldnames or [])

        if key_col == "auto":
            if "accession" in fieldnames:
                key_col_use = "accession"
            elif "pseudo_accession" in fieldnames:
                key_col_use = "pseudo_accession"
            elif fieldnames:
                key_col_use = fieldnames[0]
            else:
                key_col_use = "accession"

        for row in r:
            row2 = {kk: ("" if vv is None else str(vv)) for kk, vv in row.items()}
            primary_key = (row2.get(key_col_use) or "").strip()
            if primary_key:
                _add_meta_alias(meta, primary_key, row2)
            for extra_key in ["accession", "version", "pseudo_accession"]:
                if extra_key in row2 and row2.get(extra_key, "").strip():
                    _add_meta_alias(meta, row2[extra_key], row2)

    return meta, key_col_use, fieldnames


def local_row_to_ncbi_meta(row: Dict[str, str]) -> Dict[str, str]:
    if "accession" in row and any(
        k in row for k in ("version", "definition", "organism", "taxonomy")
    ):
        out = {k: str(row.get(k, "") or "") for k in NT_META_COLS}
        acc = (out.get("accession") or "").strip()
        if not acc:
            acc = accession_base(out.get("version") or "")
            out["accession"] = acc
        if not out.get("version") and acc:
            out["version"] = acc
        if not out.get("locus") and acc:
            out["locus"] = acc
        return out

    acc_full = (
        row.get("pseudo_accession") or row.get("accession") or row.get("version") or ""
    ).strip()
    acc_base = accession_base(acc_full)
    version  = acc_full if "." in acc_full else (row.get("version") or acc_full)

    organism = (
        row.get("scientific_name") or row.get("organism") or
        row.get("scientific_name_ascii") or ""
    ).strip()
    if not organism:
        organism = (
            row.get("species") or row.get("species_ascii") or ""
        ).strip().replace("_", " ")

    definition = (row.get("definition") or "").strip()
    if not definition:
        marker    = (row.get("marker")  or "").strip()
        voucher   = (row.get("voucher") or "").strip()
        definition = f"{marker} reference voucher:{voucher} [{organism}]".strip()

    taxonomy = (row.get("taxonomy") or "").strip()
    if not taxonomy:
        parts = []
        for k in ["kingdom", "phylum", "class", "order", "family", "genus"]:
            v = (row.get(k) or "").strip()
            if v:
                parts.append(v)
        taxonomy = "; ".join(parts)

    comment_parts = []
    for k in [
        "marker", "voucher", "station", "location", "location_ascii", "location_original",
        "sampling_date", "taxid", "bin", "public", "public_level", "notes",
    ]:
        v = (row.get(k) or "").strip()
        if v:
            comment_parts.append(f"{k}={v}")

    out = {
        "accession":          acc_base,
        "version":            version,
        "locus":              acc_base,
        "length":             (row.get("length") or row.get("seq_len") or "").strip(),
        "definition":         definition,
        "molecule_type":      (row.get("molecule_type") or "DNA").strip() or "DNA",
        "topology":           (row.get("topology") or "linear").strip() or "linear",
        "division":           (row.get("division") or "").strip(),
        "update_date":        (row.get("update_date") or "").strip(),
        "create_date":        (row.get("create_date") or row.get("sampling_date") or "").strip(),
        "organism":           organism,
        "taxonomy":           taxonomy,
        "reference_titles":   (row.get("reference_titles") or "").strip(),
        "reference_authors":  (row.get("reference_authors") or "").strip(),
        "reference_journals": (row.get("reference_journals") or "").strip(),
        "comment":            " ; ".join(comment_parts),
    }
    for k in NT_META_COLS:
        out.setdefault(k, "")
    return out


# ─────────────────────────────────────────────
# Fix-G: Entrez metadata with ENTREZ_LIMITER and proper retry
# ─────────────────────────────────────────────
def fetch_metadata_batch(
    accession_list: List[str],
    args: argparse.Namespace,
    batch_size: int = 20,
    preferred_db: str = "nucleotide",
) -> Dict[str, dict]:
    # Codex 修改: preferred_db + accession 分组 + missing 后备库重试。
    # 理由: 原代码按整批猜 nucleotide/protein；混合批次或 protein accession 不在少数前缀时
    #       会查错库，metadata 缺失。
    # 原代码摘要:
    #   db_primary = "nucleotide"
    #   for acc in batch:
    #       if acc.startswith(("NP_", "YP_", "WP_", "XP_")):
    #           db_primary = "protein"
    #   # 整批失败后再整批查 alt db
    entrez_setup(args)

    shorter_fields = [x.strip() for x in (args.shorter_metadata or "").split(",") if x.strip()]
    if not shorter_fields:
        shorter_fields = ["organism", "reference_titles", "reference_authors", "reference_journals"]

    def want(field: str) -> bool:
        return args.metadata_type == "all" or field in shorter_fields

    meta_dict: Dict[str, dict] = {}
    preferred_db = preferred_db if preferred_db in {"nucleotide", "protein"} else "nucleotide"

    def record_to_meta(record: dict) -> Optional[Tuple[str, Dict[str, str]]]:
        accession = (
            record.get("GBSeq_primary-accession") or
            record.get("GBSeq_accession-version", "").split(".")[0]
        )
        if not accession:
            return None

        meta: Dict[str, str] = {
            "accession": accession,
            "version":   record.get("GBSeq_accession-version", ""),
        }

        for field_pair in [
            ("locus",         "GBSeq_locus"),
            ("length",        "GBSeq_length"),
            ("definition",    "GBSeq_definition"),
            ("molecule_type", "GBSeq_moltype"),
            ("topology",      "GBSeq_topology"),
            ("division",      "GBSeq_division"),
            ("update_date",   "GBSeq_update-date"),
            ("create_date",   "GBSeq_create-date"),
            ("organism",      "GBSeq_organism"),
            ("taxonomy",      "GBSeq_taxonomy"),
            ("comment",       "GBSeq_comment"),
        ]:
            local_key, gb_key = field_pair
            if want(local_key):
                meta[local_key] = record.get(gb_key, "")

        ref_fields = {"reference_titles", "reference_authors", "reference_journals"}
        if args.metadata_type == "all" or any(f in shorter_fields for f in ref_fields):
            refs = record.get("GBSeq_references", []) or []
            if refs:
                titles, authors, journals = [], [], []
                for ridx, ref in enumerate(refs, start=1):
                    if want("reference_titles") and ref.get("GBReference_title"):
                        titles.append(f"Ref{ridx}: {ref['GBReference_title']}")
                    if want("reference_authors") and ref.get("GBReference_authors"):
                        authors_str = ", ".join(ref.get("GBReference_authors") or [])
                        authors.append(f"Ref{ridx}: {authors_str}")
                    if want("reference_journals") and ref.get("GBReference_journal"):
                        journals.append(f"Ref{ridx}: {ref['GBReference_journal']}")
                if titles:
                    meta["reference_titles"]   = "; ".join(titles)
                if authors:
                    meta["reference_authors"]  = "; ".join(authors)
                if journals:
                    meta["reference_journals"] = "; ".join(journals)

        return accession, meta

    def do_fetch(db_name: str, ids: List[str]) -> object:
        ENTREZ_LIMITER.acquire()
        handle = Entrez.efetch(db=db_name, id=",".join(ids), rettype="gb", retmode="xml")
        try:
            recs = Entrez.read(handle)
        finally:
            handle.close()
        return recs

    def fetch_ids_from_db(ids: List[str], db_name: str, label: str) -> Tuple[List[dict], Optional[Exception]]:
        if not ids:
            return [], None
        sleep = max(1, int(args.retry_sleep))
        last_err: Optional[Exception] = None
        for attempt in range(1, args.retries + 1):
            try:
                recs = list(do_fetch(db_name, ids))
                return recs, None
            except Exception as e:
                last_err = e
                logger.warning(
                    f"[META_FETCH] {label} db={db_name} n={len(ids)} "
                    f"attempt {attempt}/{args.retries}: {e}"
                )
                if attempt < args.retries:
                    time.sleep(sleep)
                    sleep = int(sleep * 1.6)
        return [], last_err

    grouped: Dict[str, List[str]] = {"nucleotide": [], "protein": []}
    for acc in accession_list:
        db_name = guess_entrez_db_for_accession(acc, default_db=preferred_db)
        grouped[db_name].append(acc)

    batch_counter = 0
    for db_primary in ["nucleotide", "protein"]:
        ids_for_db = grouped[db_primary]
        db_alt = "protein" if db_primary == "nucleotide" else "nucleotide"
        for i in range(0, len(ids_for_db), batch_size):
            batch_counter += 1
            batch = ids_for_db[i:i + batch_size]
            if not batch:
                continue

            recs, primary_err = fetch_ids_from_db(batch, db_primary, f"batch {batch_counter}")
            found: set = set()
            for record in recs:
                parsed = record_to_meta(record)
                if parsed is None:
                    continue
                accession, meta = parsed
                meta_dict[accession] = meta
                found.add(accession_base(accession))
                version = accession_base(meta.get("version", ""))
                if version:
                    found.add(version)

            missing_ids = [a for a in batch if accession_base(a) not in found]
            if missing_ids:
                logger.info(
                    f"[META_FETCH] trying alt db={db_alt} for {len(missing_ids)} missing "
                    f"accessions from batch {batch_counter}"
                )
                alt_recs, alt_err = fetch_ids_from_db(missing_ids, db_alt, f"batch {batch_counter} alt")
                for record in alt_recs:
                    parsed = record_to_meta(record)
                    if parsed is None:
                        continue
                    accession, meta = parsed
                    meta_dict[accession] = meta
                    found.add(accession_base(accession))
                    version = accession_base(meta.get("version", ""))
                    if version:
                        found.add(version)

                still_missing = [a for a in missing_ids if accession_base(a) not in found]
                if still_missing:
                    logger.warning(
                        f"[META_FETCH_FAIL] batch {batch_counter}: {len(still_missing)} accessions "
                        f"not found after {db_primary}/{db_alt} attempts. "
                        f"primary_err={primary_err}; alt_err={alt_err}"
                    )

            time.sleep(float(args.entrez_delay_sec))

    return meta_dict


# ─────────────────────────────────────────────
# merge with metadata
# ─────────────────────────────────────────────
def merge_rows_with_metadata(
    rows: List[Dict[str, object]],
    metadata_dict: Dict[str, dict],
    local_meta: Dict[str, Dict[str, str]],
    args: argparse.Namespace,
) -> List[Dict[str, object]]:
    if not rows:
        return rows
    out: List[Dict[str, object]] = []
    # Codex 修改: metadata length 改写入 metadata_length。
    # 理由: BLAST_HEADERS 已有 length（比对长度）；原 setdefault("length", ...) 会丢失 GenBank 序列长度。
    def add_meta_field(target: Dict[str, object], key: str, value: object) -> None:
        if key == "length" and "length" in target:
            target.setdefault("metadata_length", value)
        else:
            target.setdefault(key, value)

    for r in rows:
        rr       = dict(r)
        acc_full = parse_accession_from_sseqid(str(rr.get("sseqid", "")))
        acc      = accession_base(acc_full)

        md   = metadata_dict.get(acc)
        lrow = local_meta.get(acc_full) or local_meta.get(acc)

        if md:
            for k, v in md.items():
                if k == "accession":
                    continue
                add_meta_field(rr, k, v)

        if lrow:
            if args.local_meta_mode == "map_to_ncbi":
                lm = local_row_to_ncbi_meta(lrow)
                for k, v in lm.items():
                    if k == "accession":
                        continue
                    add_meta_field(rr, k, v)
            else:
                for k, v in lrow.items():
                    rr.setdefault(args.local_meta_prefix + k, v)

        out.append(rr)
    return out


# Codex 修改: 抽出通用分隔符写表逻辑，并新增 rows_to_csv。
# 理由: 原 rows_to_tsv 被用于写 .csv 路径，导致扩展名是 CSV 但内容是 TSV。
# 原代码保留为注释摘要:
#   def rows_to_tsv(path, rows): ... csv.DictWriter(..., delimiter="\t")
def _ordered_row_keys(
    rows: List[Dict[str, object]], default_headers: Optional[Sequence[str]] = None
) -> List[str]:
    base_headers = list(default_headers or BLAST_HEADERS)
    if not rows:
        return base_headers

    keys = list(rows[0].keys())
    ordered = [h for h in base_headers if h in keys] + [k for k in keys if k not in base_headers]
    seen = set(ordered)
    for r in rows[1:]:
        for k in r.keys():
            if k not in seen:
                ordered.append(k)
                seen.add(k)
    return ordered


def rows_to_delimited(
    path: Path,
    rows: List[Dict[str, object]],
    delimiter: str,
    default_headers: Optional[Sequence[str]] = None,
) -> None:
    ensure_dir(path.parent)
    ordered = _ordered_row_keys(rows, default_headers=default_headers)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ordered, delimiter=delimiter, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in ordered})


def rows_to_tsv(
    path: Path, rows: List[Dict[str, object]], default_headers: Optional[Sequence[str]] = None
) -> None:
    rows_to_delimited(path, rows, delimiter="\t", default_headers=default_headers)


def rows_to_csv(
    path: Path, rows: List[Dict[str, object]], default_headers: Optional[Sequence[str]] = None
) -> None:
    rows_to_delimited(path, rows, delimiter=",", default_headers=default_headers)


def debug_local_hit_examples(
    rows: List[Dict[str, object]],
    local_meta: Dict[str, Dict[str, str]],
    prefixes: List[str],
    n: int = 10,
) -> None:
    if not rows:
        return
    logger.info("[LOCAL_META_DEBUG] first parsed hit examples:")
    for i, r in enumerate(rows[:n], start=1):
        raw      = str(r.get("sseqid", ""))
        acc_full = parse_accession_from_sseqid(raw)
        acc_base = accession_base(acc_full)
        is_local = is_local_or_pseudo_accession(acc_full, local_meta, prefixes)
        logger.info(
            f"[LOCAL_META_DEBUG] {i}: sseqid={raw} | parsed={acc_full} | "
            f"base={acc_base} | is_local={is_local}"
        )



# ─────────────────────────────────────────────
# I/O discovery
# ─────────────────────────────────────────────
def find_input_files(
    input_path: str, extensions: List[str], suffix: str, mode: str
) -> List[str]:
    p = Path(input_path)
    if mode == "metadata-only":
        if p.is_file():
            return [str(p)]
        if p.is_dir():
            # Codex 修改: metadata-only 目录输入支持常见 accession-list 扩展名。
            # 理由: 原代码只查 *.txt，旧版帮助中允许 accession list，实际 .tsv/.csv 常见。
            # 原代码:
            #   return [str(x) for x in sorted(p.glob("*.txt"))]
            files: List[Path] = []
            for ext in ACCESSION_LIST_EXTENSIONS:
                files.extend(p.glob(f"*{ext}"))
            return [str(x) for x in sorted(set(files))]
        return []

    if mode == "topN-only":
        # Codex 修改: topN-only 允许直接输入已有 BLAST 表，或目录中的 BLAST 表。
        # 理由: 原逻辑按 FASTA 扩展名找输入，导致 topN-only 无法从既有 BLAST 结果生成 TopN。
        if p.is_file():
            return [str(p)]
        if p.is_dir():
            files: List[Path] = []
            for ext in list(extensions) + BLAST_TABLE_EXTENSIONS:
                files.extend(p.glob(f"*{ext}"))
            return [str(x) for x in sorted(set(files))]
        return []

    if p.is_file():
        ok = any(str(p).endswith(ext + suffix) for ext in extensions)
        return [str(p)] if ok else []

    if p.is_dir():
        files: List[str] = []
        for ext in extensions:
            files.extend(glob.glob(str(p / f"*{ext}{suffix}")))
        return sorted(files)

    return []


def setup_output_directory(output_dir: Path) -> Dict[str, Path]:
    output_dir = ensure_dir(output_dir)
    return {
        "blast":   ensure_dir(output_dir / "blast_results"),
        "topn":    ensure_dir(output_dir / "topN_results"),
        "meta":    ensure_dir(output_dir / "metadata_results"),
        "summary": ensure_dir(output_dir / "summary_results"),
        "temp":    ensure_dir(output_dir / "temp"),
    }


def _write_completion(marker: Path, args: argparse.Namespace, seq_type: Optional[str]) -> None:
    with open(marker, "w", encoding="utf-8") as f:
        f.write(f"completed_at: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"mode: {args.mode}\n")
        if seq_type:
            f.write(f"seq_type: {seq_type}\n")
        f.write(f"blast_method: {args.blast_method}\n")
        f.write(f"threads: {args.threads}\n")
        f.write(f"parallel_jobs: {args.parallel_jobs}\n")
        f.write(f"batch_size: {args.batch_size}\n")
        f.write(f"pandas: {HAVE_PANDAS}\n")
        f.write(f"entrez_tool: {args.entrez_tool}\n")
        f.write(f"api_key_set: {bool(args.api_key)}\n")


def preflight_checks(args: argparse.Namespace, files: List[str], out_dir: Path) -> None:
    issues:   List[str] = []
    warnings: List[str] = []
    # Codex 修改: 去重 fatal issue，避免多个输入文件重复报告同一缺失程序。
    # 理由: 新增按输入类型检查 BLAST+ 可执行程序后，重复信息会降低可读性。
    issue_seen = set()

    def add_issue(msg: str) -> None:
        if msg not in issue_seen:
            issue_seen.add(msg)
            issues.append(msg)

    if not files:
        add_issue("No input FASTA / input files found.")

    for fp in files[:20]:
        p = Path(fp)
        if not p.exists():
            add_issue(f"Missing input file: {p}")
        elif p.stat().st_size == 0:
            add_issue(f"Empty input file: {p}")

    if args.mode in BLAST_REQUIRED_MODES:
        # Codex 修改: BLAST 预检按实际序列类型和 online_engine 检查 blastn/blastp。
        # 理由: 原代码只在 local 模式固定检查 nucleotide/blastn；默认 blast_remote 缺 blastn
        #       会到运行阶段才失败，并且旧 main 还返回 0。
        # 原代码:
        #   if args.blast_method == "local":
        #       program = args.blast_exe if args.blast_exe != "auto" else get_blast_program("nucleotide")
        #       ...
        for fp in files[:20]:
            seq_type = detect_sequence_type(fp)
            program = args.blast_exe if args.blast_exe != "auto" else get_blast_program(seq_type)
            needs_blast_plus = args.blast_method == "local" or (
                args.blast_method == "online" and args.online_engine == "blast_remote"
            )
            if needs_blast_plus and not which(program):
                add_issue(f"BLAST+ executable not found in PATH for {seq_type} input: {program}")

            if args.blast_method == "local":
                db = get_blast_db(seq_type, args)
                db_p = Path(str(db))
                suffixes = (
                    [".pin", ".phr", ".psq", ".00.pin", ".00.phr", ".00.psq"]
                    if seq_type == "protein"
                    else [".nin", ".nhr", ".nsq", ".00.nin", ".00.nhr", ".00.nsq"]
                )
                if db not in {"nt", "nr"} and not db_p.exists() and not any(
                    Path(str(db) + s).exists() for s in suffixes
                ):
                    warnings.append(f"Local BLAST DB path could not be verified directly: {db}")

    if args.local_meta_tsv:
        lp = Path(args.local_meta_tsv)
        if not lp.exists():
            add_issue(f"--local_meta_tsv not found: {lp}")
        elif lp.stat().st_size == 0:
            add_issue(f"--local_meta_tsv is empty: {lp}")

    if args.mode in {"blast2meta", "blast2topN2metadata", "metadata-only"}:
        # Codex 修改: 对 placeholder email 给出预警。
        # 理由: Entrez 请求需要真实 email，默认占位符容易造成合规和限流问题。
        if not args.no_ncbi_metadata and args.email in {"", "you@example.com"}:
            warnings.append(
                "NCBI metadata fetch is enabled but --email is empty or the placeholder you@example.com."
            )

    if args.use_snakemake and not which("snakemake"):
        warnings.append("snakemake not found; workflow will fall back to plain Python mode.")

    out_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "fatal_issues":  issues,
        "warnings":      warnings,
        "n_input_files": len(files),
        "input_preview": files[:20],
        "output_dir":    str(out_dir),
        "timestamp":     time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(out_dir / "preflight_report.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)

    for w in warnings:
        logger.warning(f"[PREFLIGHT] {w}")
    for e in issues:
        logger.error(f"[PREFLIGHT] {e}")

    if issues:
        raise SystemExit("Preflight checks failed. See preflight_report.json")


# ─────────────────────────────────────────────
# main per-FASTA (Fix-H: csv_output always written)
# ─────────────────────────────────────────────
def process_one_fasta(
    fasta_file: str,
    args: argparse.Namespace,
    dirs: Dict[str, Path],
) -> bool:
    fasta_path = Path(fasta_file)
    # Codex 修改: process_one_fasta 返回 bool，并识别 topN-only 的 BLAST 表输入。
    # 理由: 原代码返回 None，失败路径只 log 后 return，主程序无法产生非零退出码；
    #       topN-only 也无法直接处理已有 BLAST TSV。
    # 原签名:
    #   def process_one_fasta(... ) -> None:
    fasta_exts = [x.strip() for x in args.extension_fasta.split(",") if x.strip()]
    input_is_blast_table = args.mode == "topN-only" and not looks_like_fasta_path(fasta_path, fasta_exts)
    file_stem = (
        args.out_stem
        if args.out_stem
        else strip_blast_suffix(fasta_path.stem) if input_is_blast_table else fasta_path.stem
    )

    blast_output_merged = dirs["blast"]   / f"{file_stem}_blast.tsv"
    topn_output         = dirs["topn"]    / f"{file_stem}_topN.tsv"
    metadata_output     = dirs["meta"]    / f"{file_stem}_metadata.tsv"
    summary_output      = dirs["summary"] / f"{file_stem}_summary.tsv"
    csv_output          = Path(args.output) / f"{file_stem}_results.csv"
    completion_marker   = Path(args.output) / f"{file_stem}_completed.txt"

    if args.resume and is_file_complete(csv_output) and completion_marker.exists():
        logger.info(f"[SKIP] {fasta_file} (resume: outputs exist)")
        return True

    blast_needed = args.mode in ["blast-only", "blast2topN", "blast2meta", "blast2topN2metadata"]
    topn_needed  = args.mode in ["topN-only", "blast2topN", "blast2topN2metadata"]
    meta_needed  = args.mode in ["blast2meta", "blast2topN2metadata", "metadata-only"]

    seq_type = None
    # Codex 修改: 只在真正要跑 BLAST 时检测序列类型。
    # 理由: 原代码在 topN-only 也尝试把 BLAST TSV 当 FASTA 检测，容易误导日志。
    # 原代码:
    #   if args.mode != "metadata-only":
    if blast_needed:
        seq_type = detect_sequence_type(fasta_path)
        logger.info(f"[SEQTYPE] {fasta_path.name}: {seq_type}")

    local_meta: Dict[str, Dict[str, str]] = {}
    if args.local_meta_tsv:
        local_meta, local_meta_key_used, fieldnames = load_local_meta(
            args.local_meta_tsv, args.local_meta_key
        )
        logger.info(
            f"[LOCAL_META] loaded rows={len(local_meta)} "
            f"key_requested={args.local_meta_key} key_used={local_meta_key_used}"
        )
        if fieldnames:
            logger.info(f"[LOCAL_META] columns={','.join(fieldnames)}")
        if len(local_meta) == 0:
            logger.warning(
                "[LOCAL_META] Loaded 0 rows. Check --local_meta_tsv, and whether the key "
                "matches the file. For local_nt_metadata.tsv the best key is usually "
                "accession or auto."
            )

    pseudo_prefixes = [
        x.strip() for x in (args.pseudo_acc_prefixes or "").split(",") if x.strip()
    ]

    # Fix-2/3: derive col_names once so parser and writers stay consistent
    _, col_names = build_outfmt_string(args.outfmt_blast)
    # Codex 修改: qblast 后端禁止自定义列顺序。
    # 理由: Biopython qblast 解析 XML 后固定写 14 列，不能遵守 BLAST+ outfmt 的自定义列顺序。
    if (
        blast_needed
        and args.blast_method == "online"
        and args.online_engine == "biopython"
        and col_names != BLAST_HEADERS
    ):
        logger.error(
            "[QBLAST] --online_engine biopython writes the built-in 14 BLAST columns; "
            "custom --outfmt_blast column order is only supported with BLAST+ engines."
        )
        return False

    # ── 1) BLAST ──────────────────────────────
    if blast_needed:
        if args.resume and is_file_complete(blast_output_merged):
            logger.info(f"[BLAST] merged exists, skip: {blast_output_merged.name}")
        else:
            batch_size = int(args.batch_size) if int(args.batch_size) > 0 else 50
            batches    = iter_fasta_batches(fasta_path, dirs["temp"] / file_stem, batch_size)
            if not batches:
                logger.warning(f"[EMPTY_FASTA] {fasta_file}")
                return False

            blast_files: List[Path] = []
            if args.blast_method == "local":
                parallel_jobs, threads_per_job = compute_parallel_params(args, mode="local")
                logger.info(
                    f"[LOCAL_SCHED] total_threads={args.threads} "
                    f"parallel_jobs={parallel_jobs} threads_per_job={threads_per_job}"
                )
                def runner(b: BatchInfo) -> Tuple[str, bool, Path, List[str]]:
                    return run_local_blast_one_batch(
                        b, seq_type, args, dirs["blast"] / file_stem, threads_per_job, col_names
                    )
            else:
                parallel_jobs, _ = compute_parallel_params(args, mode="online")
                logger.info(
                    f"[ONLINE_SCHED] parallel_jobs={parallel_jobs} "
                    f"engine={args.online_engine} batch_size={args.batch_size} "
                    f"delay={args.online_delay_sec}s"
                )
                if args.online_engine == "blast_remote":
                    def runner(b: BatchInfo) -> Tuple[str, bool, Path, List[str]]:
                        return run_remote_blast_one_batch(
                            b, seq_type, args, dirs["blast"] / file_stem, col_names
                        )
                else:
                    def runner(b: BatchInfo) -> Tuple[str, bool, Path, List[str]]:
                        return run_biopython_qblast_one_batch(
                            b, seq_type, args, dirs["blast"] / file_stem, col_names
                        )

            with cf.ThreadPoolExecutor(
                max_workers=min(parallel_jobs, len(batches))
            ) as ex:
                futs = [ex.submit(runner, b) for b in batches]
                for fut in cf.as_completed(futs):
                    bid, ok, outp, _ = fut.result()
                    if ok:
                        blast_files.append(outp)
                    else:
                        logger.error(f"[BLAST_FAIL] batch={bid}")

            if not blast_files:
                logger.error(f"[BLAST] no successful batches for {fasta_file}")
                return False

            merge_blast_tsv_text(blast_files, blast_output_merged)
            logger.info(f"[BLAST] merged -> {blast_output_merged}")

        # Fix-H.1: blast-only mode must write csv_output
        if args.mode == "blast-only":
            merged_rows = parse_blast_results_py(blast_output_merged, col_names)
            if HAVE_PANDAS and pd is not None:
                rows_to_df(merged_rows, default_headers=col_names).to_csv(csv_output, index=False)
            else:
                rows_to_csv(csv_output, merged_rows, default_headers=col_names)
            logger.info(f"[CSV] wrote -> {csv_output}")
            _write_completion(completion_marker, args, seq_type)
            return True

    merged_rows: List[Dict[str, object]] = []
    if blast_needed and blast_output_merged.exists():
        merged_rows = parse_blast_results_py(blast_output_merged, col_names)
    elif args.mode == "topN-only":
        # Codex 修改: topN-only 从已有 BLAST 表加载 merged_rows。
        # 理由: 原代码中 blast_needed=False，merged_rows 永远为空，导致静默输出空 TopN。
        # 原代码:
        #   merged_rows = []
        #   if blast_needed and blast_output_merged.exists(): ...
        blast_source = fasta_path if input_is_blast_table else blast_output_merged
        if not is_file_complete(blast_source):
            logger.error(
                f"[TOPN] BLAST results not found for topN-only mode: {blast_source}. "
                "Use a BLAST TSV as -i or keep the prior output under blast_results/{stem}_blast.tsv."
            )
            return False
        merged_rows = parse_blast_results_py(blast_source, col_names)
        logger.info(f"[TOPN] loaded BLAST table -> {blast_source}")

    if local_meta and merged_rows:
        debug_local_hit_examples(merged_rows, local_meta, pseudo_prefixes, n=8)

    # ── 2) TopN ──────────────────────────────
    top_rows: Optional[List[Dict[str, object]]] = None
    if topn_needed:
        if args.resume and is_file_complete(topn_output):
            logger.info(f"[TOPN] exists, skip: {topn_output.name}")
            top_rows = parse_blast_results_py(topn_output, col_names)
        else:
            if not merged_rows:
                logger.warning(f"[TOPN] no hits for {fasta_file}")
                top_rows = []
            else:
                top_rows = get_topN_hits_py(merged_rows, int(args.topN), args.sort_by)
            rows_to_tsv(topn_output, top_rows or [], default_headers=col_names)
            logger.info(f"[TOPN] wrote -> {topn_output}")

        # Fix-H.2: topN-only mode must write csv_output directly
        if args.mode == "topN-only":
            if HAVE_PANDAS and pd is not None:
                rows_to_df(top_rows or [], default_headers=col_names).to_csv(csv_output, index=False)
            else:
                rows_to_csv(csv_output, top_rows or [], default_headers=col_names)
            logger.info(f"[CSV] wrote -> {csv_output}")
            _write_completion(completion_marker, args, seq_type)
            return True

        if args.mode == "blast2topN":
            # Codex 修改: blast2topN 也写最终 *_results.csv。
            # 理由: 原代码只写 topN_results/*.tsv，Snakemake rule all 和 --resume 期待 CSV 哨兵。
            if HAVE_PANDAS and pd is not None:
                rows_to_df(top_rows or [], default_headers=col_names).to_csv(csv_output, index=False)
            else:
                rows_to_csv(csv_output, top_rows or [], default_headers=col_names)
            logger.info(f"[CSV] wrote -> {csv_output}")
            _write_completion(completion_marker, args, seq_type)
            return True

    # ── 3) Metadata fetch + merge ─────────────
    if meta_needed:
        metadata_dict: Dict[str, dict] = {}
        # Codex 修改: metadata-only 读取 accession list；其他模式仍从 BLAST/TopN 行提取。
        # 理由: 原代码在 metadata-only 下直接 accessions=[]，该模式实际空跑。
        # 原代码:
        #   if args.mode == "metadata-only":
        #       accessions: List[str] = []
        if args.mode == "metadata-only":
            accessions = read_accession_list(fasta_path)
            logger.info(f"[META_ONLY] loaded accessions={len(accessions)} from {fasta_path}")
        elif args.mode == "blast2meta":
            accessions = extract_accession_numbers_from_rows(merged_rows)
        else:
            accessions = extract_accession_numbers_from_rows(top_rows or [])

        if args.mode == "metadata-only" and local_meta:
            # Codex 修改: metadata-only + local_meta_tsv 时可直接输出本地 metadata。
            # 理由: 旧逻辑只有 BLAST summary merge 才能使用 local_meta，metadata-only 无法替代本地库用法。
            for acc in accessions:
                acc_full = (acc or "").strip()
                acc_base = accession_base(acc_full)
                lrow = local_meta.get(acc_full) or local_meta.get(acc_base)
                if not lrow:
                    continue
                if args.local_meta_mode == "map_to_ncbi":
                    lm = local_row_to_ncbi_meta(lrow)
                    key = lm.get("accession") or acc_base
                    metadata_dict[key] = lm
                else:
                    metadata_dict[acc_base] = {
                        "accession": acc_base,
                        **{args.local_meta_prefix + k: v for k, v in lrow.items()},
                    }

        if not args.no_ncbi_metadata:
            accessions_real = [
                a for a in accessions
                if not is_local_or_pseudo_accession(a, local_meta, pseudo_prefixes)
            ]
            logger.info(
                f"[META] accessions(total)={len(accessions)} "
                f"real_for_entrez={len(accessions_real)} "
                f"local_or_pseudo_skipped={len(accessions) - len(accessions_real)}"
            )

            if args.resume and is_file_complete(metadata_output):
                logger.info(f"[META] exists, load: {metadata_output.name}")
                with open(metadata_output, "r", encoding="utf-8", errors="ignore") as f:
                    reader = csv.DictReader(f, delimiter="\t")
                    for r in reader:
                        a = str(r.get("accession", "")).split(".")[0]
                        if a:
                            metadata_dict[a] = dict(r)

            missing = [a for a in accessions_real if a not in metadata_dict]
            if missing:
                preferred_db = "protein" if seq_type == "protein" else "nucleotide"
                metadata_new = fetch_metadata_batch(
                    missing, args, batch_size=20, preferred_db=preferred_db
                )
                metadata_dict.update(metadata_new)

            if metadata_dict:
                rows_to_tsv(metadata_output, list(metadata_dict.values()), default_headers=NT_META_COLS)
                logger.info(f"[META] wrote -> {metadata_output}")
        else:
            logger.info("[META] --no_ncbi_metadata ON, skip Entrez efetch.")
            if metadata_dict:
                rows_to_tsv(metadata_output, list(metadata_dict.values()), default_headers=NT_META_COLS)
                logger.info(f"[META] wrote local metadata -> {metadata_output}")

        # Fix-H.3: metadata-only mode must write csv_output directly
        if args.mode == "metadata-only":
            if HAVE_PANDAS and pd is not None:
                rows_to_df(list(metadata_dict.values()), default_headers=NT_META_COLS).to_csv(csv_output, index=False)
            else:
                rows_to_csv(csv_output, list(metadata_dict.values()), default_headers=NT_META_COLS)
            logger.info(f"[CSV] wrote -> {csv_output}")
            _write_completion(completion_marker, args, seq_type)
            return True

        base_rows = merged_rows if args.mode == "blast2meta" else (top_rows or [])
        merged    = merge_rows_with_metadata(base_rows, metadata_dict, local_meta, args)

        rows_to_tsv(summary_output, merged)
        # Fix-H.4: blast2*/no-pandas fallback must write csv_output directly
        if HAVE_PANDAS and pd is not None:
            rows_to_df(merged).to_csv(csv_output, index=False)
        else:
            rows_to_csv(csv_output, merged)

        logger.info(f"[SUMMARY] wrote -> {summary_output}")
        logger.info(f"[CSV] wrote -> {csv_output}")

    _write_completion(completion_marker, args, seq_type)
    return True



# ─────────────────────────────────────────────
# Fix-1 / Fix-5: Snakemake integration
# ─────────────────────────────────────────────
def _write_text(path: Path, text: str) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def _write_samples_tsv(path: Path, files: List[str]) -> None:
    ensure_dir(path.parent)
    counts: Dict[str, int] = {}
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sample", "fasta"], delimiter="\t")
        w.writeheader()
        for fp in files:
            stem          = Path(fp).stem
            counts[stem]  = counts.get(stem, 0) + 1
            sample        = stem if counts[stem] == 1 else f"{stem}_{counts[stem]}"
            w.writerow({"sample": sample, "fasta": str(Path(fp).resolve())})


def _snakefile_text() -> str:
    """
    Fix-5: child jobs now receive --args_json instead of a long manually-assembled
    CLI.  Per-job overrides (--input, --output, --out_stem, --resume) are still
    passed explicitly so Snakemake wildcards work correctly.
    """
    return r'''
import os, csv, shlex

OUTDIR        = config["outdir"]
PY            = config.get("python_exe", "python")
SCRIPT        = config["script"]
ARGS_JSON     = config["args_json"]
LATENCY_WAIT  = int(config.get("latency_wait", 60))
THREADS_PER_JOB  = int(config["threads_per_job"])
MEM_MB_PER_JOB   = int(config.get("mem_mb_per_job", 0))
JOB_SLOTS_PER_JOB = int(config.get("job_slots_per_job", 1))


def load_samples(tsv):
    rows = []
    with open(tsv, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rows.append(row)
    return rows


SAMPLES    = load_samples(config["samples_tsv"])
BY_SAMPLE  = {r["sample"]: r for r in SAMPLES}
SAMPLE_IDS = [r["sample"] for r in SAMPLES]


rule all:
    input:
        expand(os.path.join(OUTDIR, "{sample}", "{sample}_results.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(OUTDIR, "{sample}", "{sample}_completed.txt"), sample=SAMPLE_IDS)


rule run_one:
    input:
        lambda wc: BY_SAMPLE[wc.sample]["fasta"]
    output:
        csv  = os.path.join(OUTDIR, "{sample}", "{sample}_results.csv"),
        mark = os.path.join(OUTDIR, "{sample}", "{sample}_completed.txt")
    threads:
        THREADS_PER_JOB
    resources:
        mem_mb    = MEM_MB_PER_JOB,
        job_slots = JOB_SLOTS_PER_JOB
    params:
        outdir       = OUTDIR,
        sample_outdir = lambda wc: os.path.join(OUTDIR, wc.sample),
        py           = PY,
        script       = SCRIPT,
        args_json    = ARGS_JSON,
        latency_wait = LATENCY_WAIT
    log:
        os.path.join(OUTDIR, "{sample}", "snakemake.run.log")
    shell:
        r"""
        set -euo pipefail
        # Codex 修改: Snakemake shell 参数使用 :q 引用，并把 sample_outdir 做成参数。
        # 理由: 原代码未引用路径，输入/输出/脚本路径包含空格或特殊字符时会失败。
        # 原代码:
        #   mkdir -p {params.outdir}/{wildcards.sample}
        #   {params.py} {params.script} -i {input} -o {params.outdir}/{wildcards.sample} ...
        mkdir -p {params.sample_outdir:q}
        {params.py:q} {params.script:q} \
            -i {input:q} \
            -o {params.sample_outdir:q} \
            --args_json {params.args_json:q} \
            --out_stem {wildcards.sample:q} \
            --resume \
            --no_smk \
            > {log:q} 2>&1
        test -s {output.csv:q}
        test -s {output.mark:q}
        """
'''


def _snakemake_error_indicates_lock(stderr: str, stdout: str) -> bool:
    txt     = f"{stderr}\n{stdout}".lower()
    needles = [
        "lockexception", "directory cannot be locked",
        "error: directory cannot be locked",
        "unable to acquire lock", "failed to lock",
    ]
    return any(x in txt for x in needles)


def _run_snakemake_unlock(
    smk: str, snakefile: Path, config_json: Path, wf_dir: Path
) -> Tuple[int, str, str]:
    unlock_cmd = [smk, "-s", str(snakefile), "--configfile", str(config_json), "--unlock"]
    logger.warning("[SMK] Detected lock-related failure, trying automatic --unlock once.")
    return run_cmd(unlock_cmd, cwd=wf_dir, retries=1, retry_sleep=1)


@dataclass
class SnakemakePlan:
    smk_cores:              int
    smk_jobs:               int
    threads_per_job:        int
    total_mem_mb:           int
    mem_mb_per_job:         int
    total_job_slots:        int
    job_slots_per_job:      int
    detected_cpu_count:     int
    detected_mem_total_mb:  int
    detected_mem_available_mb: int
    detected_cpu_percent:   float
    detected_loadavg1:      float
    detected_load_ratio_1m: float
    detected_mem_pressure:  float
    busy_score:             float
    busy_guard_applied:     bool
    busy_cpu_cap:           int
    busy_mem_cap_mb:        int


def _derive_busy_caps(
    host: HostState,
    requested_cores: int,
    requested_mem_mb: int,
    args: argparse.Namespace,
) -> Tuple[int, int, bool]:
    if not args.smk_busy_guard:
        return requested_cores, requested_mem_mb, False

    score = host.busy_score
    if score >= 1.10:
        cpu_frac, avail_mem_frac = 0.15, 0.35
    elif score >= 0.90:
        cpu_frac, avail_mem_frac = 0.25, 0.45
    elif score >= 0.75:
        cpu_frac, avail_mem_frac = 0.35, 0.55
    elif score >= 0.60:
        cpu_frac, avail_mem_frac = 0.50, 0.65
    elif score >= 0.45:
        cpu_frac, avail_mem_frac = 0.65, 0.75
    else:
        cpu_frac      = min(0.90, max(0.30, float(args.smk_idle_cpu_fraction)))
        avail_mem_frac = min(0.90, max(0.35, float(args.smk_idle_avail_mem_fraction)))

    user_cap_avail = float(args.smk_max_available_mem_fraction)
    user_cap_avail = min(0.95, max(0.10, user_cap_avail))
    avail_mem_frac = min(avail_mem_frac, user_cap_avail)

    min_cores  = max(1, safe_int(args.smk_busy_guard_min_cores,  1))
    min_mem_mb = max(1024, safe_int(args.smk_busy_guard_min_mem_mb, 4096))

    cpu_cap    = max(min_cores,  int(host.cpu_count      * cpu_frac))
    mem_cap_mb = max(min_mem_mb, int(host.mem_available_mb * avail_mem_frac))

    final_cores  = max(min_cores,  min(requested_cores,  cpu_cap))
    final_mem_mb = max(min_mem_mb, min(requested_mem_mb, mem_cap_mb))
    busy_applied = final_cores < requested_cores or final_mem_mb < requested_mem_mb
    return final_cores, final_mem_mb, busy_applied


def build_snakemake_plan(args: argparse.Namespace, n_files: int) -> SnakemakePlan:
    host = detect_host_state()
    cpu_count = host.cpu_count

    requested_total_threads = resolve_total_threads(args)
    requested_smk_cores     = safe_int(args.smk_cores, 0)
    smk_cores = (
        max(1, min(requested_smk_cores, cpu_count))
        if requested_smk_cores > 0 else requested_total_threads
    )
    smk_cores = max(1, min(smk_cores, cpu_count))

    total_mem_mb = safe_int(args.smk_total_mem_mb, 0)
    if total_mem_mb <= 0:
        if args.smk_auto_mem:
            base_total = int(
                host.mem_total_mb * min(0.98, max(0.05, float(args.smk_mem_fraction)))
            )
            base_avail = int(
                host.mem_available_mb * min(0.95, max(0.05, float(args.smk_max_available_mem_fraction)))
            )
            total_mem_mb = max(1024, min(base_total, base_avail))
        else:
            total_mem_mb = max(1024, int(host.mem_available_mb * 0.8))

    busy_cpu_cap  = smk_cores
    busy_mem_cap_mb = total_mem_mb
    busy_applied  = False
    if args.smk_busy_guard:
        smk_cores, total_mem_mb, busy_applied = _derive_busy_caps(host, smk_cores, total_mem_mb, args)
        busy_cpu_cap    = smk_cores
        busy_mem_cap_mb = total_mem_mb

    user_jobs           = safe_int(args.smk_jobs, 0)
    user_threads_per_job = safe_int(args.smk_threads_per_job, 0)
    user_mem_per_job    = safe_int(args.smk_mem_mb_per_job, 0)

    auto_min_threads = max(1,    safe_int(args.smk_auto_min_threads_per_job, 4))
    auto_min_mem_mb  = max(1024, safe_int(args.smk_auto_min_mem_mb_per_job, 24576))
    auto_max_jobs    = max(1,    safe_int(args.smk_auto_max_jobs, 8))

    if user_jobs > 0:
        smk_jobs = user_jobs
    else:
        jobs_by_cpu = max(1, smk_cores // auto_min_threads)
        jobs_by_mem = max(1, total_mem_mb // auto_min_mem_mb) if total_mem_mb > 0 else jobs_by_cpu
        smk_jobs    = max(1, min(max(1, n_files), jobs_by_cpu, jobs_by_mem, auto_max_jobs))

    smk_jobs = max(1, min(smk_jobs, max(1, n_files)))

    if user_threads_per_job > 0:
        threads_per_job = user_threads_per_job
    else:
        threads_per_job = max(1, smk_cores // smk_jobs)
        if threads_per_job < auto_min_threads and smk_cores >= auto_min_threads:
            threads_per_job = auto_min_threads

    threads_per_job = max(1, min(threads_per_job, smk_cores))

    max_jobs_by_threads = max(1, smk_cores // threads_per_job)
    if smk_jobs > max_jobs_by_threads:
        smk_jobs = max_jobs_by_threads

    if user_mem_per_job > 0:
        mem_mb_per_job = user_mem_per_job
    elif total_mem_mb > 0:
        mem_mb_per_job = max(1024, int(total_mem_mb / smk_jobs))
    else:
        mem_mb_per_job = 0

    if total_mem_mb > 0 and mem_mb_per_job > 0:
        max_jobs_by_mem = max(1, total_mem_mb // mem_mb_per_job)
        if smk_jobs > max_jobs_by_mem:
            smk_jobs = max_jobs_by_mem
            if user_threads_per_job <= 0:
                threads_per_job = max(1, min(smk_cores, smk_cores // smk_jobs))
            if user_mem_per_job <= 0:
                mem_mb_per_job = max(1024, int(total_mem_mb / smk_jobs))

    total_job_slots   = safe_int(args.smk_total_job_slots, 0)
    if total_job_slots <= 0:
        total_job_slots = smk_jobs
    job_slots_per_job = max(1, safe_int(args.smk_job_slots_per_job, 1))

    return SnakemakePlan(
        smk_cores=smk_cores,
        smk_jobs=smk_jobs,
        threads_per_job=threads_per_job,
        total_mem_mb=total_mem_mb,
        mem_mb_per_job=mem_mb_per_job,
        total_job_slots=total_job_slots,
        job_slots_per_job=job_slots_per_job,
        detected_cpu_count=host.cpu_count,
        detected_mem_total_mb=host.mem_total_mb,
        detected_mem_available_mb=host.mem_available_mb,
        detected_cpu_percent=host.cpu_percent,
        detected_loadavg1=host.loadavg1,
        detected_load_ratio_1m=host.load_ratio_1m,
        detected_mem_pressure=host.mem_pressure,
        busy_score=host.busy_score,
        busy_guard_applied=busy_applied,
        busy_cpu_cap=busy_cpu_cap,
        busy_mem_cap_mb=busy_mem_cap_mb,
    )


def detect_snakemake_cli_capabilities(smk: str) -> Dict[str, object]:
    caps: Dict[str, object] = {"retry_flag": None, "help_ok": False, "version": ""}

    rc, so, se = run_cmd([smk, "--help"], retries=1, retry_sleep=1)
    help_txt = (so or "") + "\n" + (se or "")
    if rc == 0 and help_txt:
        caps["help_ok"] = True
        if "--retries" in help_txt:
            caps["retry_flag"] = "--retries"
        elif "--restart-times" in help_txt:
            caps["retry_flag"] = "--restart-times"
        elif re.search(r"(^|\s)-T([,\s]|$)", help_txt):
            caps["retry_flag"] = "-T"

    vrc, vso, vse = run_cmd([smk, "--version"], retries=1, retry_sleep=1)
    if vrc == 0:
        caps["version"] = (vso or vse or "").strip()
    return caps


def write_smk_plan_json(
    path: Path, plan: SnakemakePlan, caps: Dict[str, object]
) -> None:
    payload = {
        "snakemake_plan": asdict(plan),
        "snakemake_cli":  caps,
        "timestamp":      time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    _write_text(path, json.dumps(payload, indent=2, ensure_ascii=False))


def _args_to_json_safe(args: argparse.Namespace) -> Dict:
    """
    Convert args namespace to a JSON-serialisable dict.
    Non-serialisable values are coerced to str as a fallback.
    """
    d = {}
    for k, v in vars(args).items():
        try:
            json.dumps(v)
            d[k] = v
        except (TypeError, ValueError):
            d[k] = str(v)
    return d


def run_with_snakemake(
    args: argparse.Namespace, files: List[str], out_dir: Path
) -> int:
    smk = which("snakemake")
    if not smk:
        # Codex 修改: --smk_dry_run 在缺少 snakemake 时不再落回真实执行。
        # 理由: 原 fallback 会忽略 dry-run 语义，可能意外启动 BLAST/metadata 请求。
        if args.smk_dry_run:
            logger.error("[SMK] --smk_dry_run requested but snakemake is not installed; not running fallback workflow.")
            return 1
        logger.warning("[SMK] snakemake not found, fallback to python workflow.")
        dirs = setup_output_directory(out_dir)
        # Codex 修改: fallback Python workflow 汇总失败数并返回非零。
        # 理由: 原代码 catch 后仍 return 0，调度系统无法感知失败。
        failed = 0
        for fp in files:
            try:
                if not process_one_fasta(fp, args, dirs):
                    failed += 1
            except Exception as e:
                failed += 1
                logger.error(f"[FILE_FAIL] {fp}: {e}", exc_info=True)
        return 1 if failed else 0

    wf_dir      = ensure_dir(out_dir / "_snakemake_workflow")
    samples_tsv = wf_dir / "samples.tsv"
    config_json = wf_dir / "config.json"
    snakefile   = wf_dir / "Snakefile"
    plan_json   = wf_dir / "snakemake_plan.json"

    # Fix-1: serialise FULL args namespace so child jobs auto-inherit everything
    child_args_json = wf_dir / "child_args.json"
    child_args_dict = _args_to_json_safe(args)
    # Override fields that must differ in child processes
    child_args_dict["use_snakemake"] = False
    # threads_per_job will be overridden below once plan is built; placeholder for now
    _write_text(child_args_json, json.dumps(child_args_dict, indent=2, ensure_ascii=False))

    _write_samples_tsv(samples_tsv, files)
    _write_text(snakefile, _snakefile_text())

    plan = build_snakemake_plan(args, len(files))
    caps = detect_snakemake_cli_capabilities(smk)

    # Update the child args JSON with the actual per-job thread count decided by plan.
    # Fix-A: also inject active-job count so each child process can compute its
    # fair share of the global NCBI rate budget (3 or 10 req/sec ÷ smk_jobs).
    child_args_dict["threads"] = int(plan.threads_per_job)
    child_args_dict["_smk_active_jobs"] = int(plan.smk_jobs)
    _write_text(child_args_json, json.dumps(child_args_dict, indent=2, ensure_ascii=False))

    write_smk_plan_json(plan_json, plan, caps)

    logger.info(
        "[SMK_PLAN] detected_cpu=%s total_mem_mb=%s avail_mem_mb=%s "
        "cpu_percent=%.1f load1=%.2f load_ratio=%.2f mem_pressure=%.2f busy_score=%.2f "
        "busy_guard=%s cores=%s jobs=%s threads_per_job=%s "
        "total_mem_mb=%s mem_mb_per_job=%s total_job_slots=%s job_slots_per_job=%s",
        plan.detected_cpu_count, plan.detected_mem_total_mb, plan.detected_mem_available_mb,
        plan.detected_cpu_percent, plan.detected_loadavg1,
        plan.detected_load_ratio_1m, plan.detected_mem_pressure, plan.busy_score,
        plan.busy_guard_applied, plan.smk_cores, plan.smk_jobs, plan.threads_per_job,
        plan.total_mem_mb, plan.mem_mb_per_job, plan.total_job_slots, plan.job_slots_per_job,
    )
    logger.info("[SMK_CLI] version=%s retry_flag=%s", caps.get("version", ""), caps.get("retry_flag"))

    # Fix-5: config now carries args_json path instead of a long base_args list
    cfg: Dict[str, object] = {
        "outdir":           str(out_dir.resolve()),
        "python_exe":       sys.executable,
        "script":           str(Path(__file__).resolve()),
        "samples_tsv":      str(samples_tsv.resolve()),
        "args_json":        str(child_args_json.resolve()),
        "threads_per_job":  int(plan.threads_per_job),
        "mem_mb_per_job":   int(plan.mem_mb_per_job),
        "job_slots_per_job": int(plan.job_slots_per_job),
        "latency_wait":     int(args.smk_latency_wait),
    }
    _write_text(config_json, json.dumps(cfg, indent=2, ensure_ascii=False))

    # Fix-B: In Snakemake v8+ `--cores` / `-c` and `--jobs` / `-j` are strict
    # aliases.  Passing both causes: "error: not allowed with argument --jobs/-j".
    # We only pass `--cores`; concurrency is controlled via --resources job_slots.
    cmd = [
        smk,
        "-s", str(snakefile),
        "--configfile", str(config_json),
        "--cores", str(plan.smk_cores),
        "--latency-wait", str(int(args.smk_latency_wait)),
    ]

    retry_flag = str(caps.get("retry_flag") or "").strip()
    if retry_flag:
        cmd += [retry_flag, str(int(args.smk_restart_times))]
    else:
        logger.warning(
            "[SMK] Could not detect a supported retry flag from snakemake --help; "
            "skipping CLI retry option."
        )

    if args.smk_rerun_incomplete:
        cmd.append("--rerun-incomplete")
    if args.smk_keep_going:
        cmd.append("--keep-going")
    if args.smk_printshellcmds:
        cmd.append("--printshellcmds")
    if args.smk_dry_run:
        cmd.append("--dry-run")

    resources_chunks: List[str] = []
    if plan.total_mem_mb > 0 and plan.mem_mb_per_job > 0:
        resources_chunks.append(f"mem_mb={plan.total_mem_mb}")
    if plan.total_job_slots > 0 and plan.job_slots_per_job > 0:
        resources_chunks.append(f"job_slots={plan.total_job_slots}")
    if resources_chunks:
        cmd += ["--resources"] + resources_chunks

    if args.smk_profile:
        cmd += ["--profile", args.smk_profile]
    if args.smk_extra:
        cmd += _shlex_split(args.smk_extra)

    rc, so, se = run_cmd(cmd, cwd=wf_dir, retries=1, retry_sleep=1)
    if rc != 0 and args.smk_unlock_on_lock_exception and _snakemake_error_indicates_lock(se, so):
        urc, uso, use = _run_snakemake_unlock(smk, snakefile, config_json, wf_dir)
        if urc == 0:
            logger.info("[SMK] unlock success, rerunning workflow.")
            rc, so, se = run_cmd(cmd, cwd=wf_dir, retries=1, retry_sleep=1)
        else:
            logger.error(f"[SMK_UNLOCK_FAIL] rc={urc}\n{use[-2000:]}")
            return urc

    if rc != 0:
        logger.error(
            f"[SMK_FAIL] rc={rc}\nSTDOUT_tail:\n{so[-3000:]}\nSTDERR_tail:\n{se[-3000:]}"
        )
        return rc

    logger.info("[SMK] Completed.")
    return 0



# ─────────────────────────────────────────────
# CLI (Fix-J: updated epilog examples)
# ─────────────────────────────────────────────
# Codex 修改(help-short): 新增短帮助文本生成函数，供 -h 和无参数错误提示复用。
# 理由: 日常运行只需要核心用法；完整高级参数仍由 argparse 的 --help 输出。
# 原代码摘要: 原来 -h/--help 都直接使用 argparse 自动生成的完整帮助。
def short_help_text(script_name: str) -> str:
    return f"""usage:
  python {script_name} -i INPUT -o OUT [options]

常用一键运行:
  # 默认 remote / online BLAST
  python {script_name} -i query.fasta -o OUT

  # 本地 BLAST 数据库
  python {script_name} -i query.fasta -o OUT -bm local --blast_db_nt my_db

  # 不使用 Snakemake，直接普通 Python 运行
  python {script_name} -i query.fasta -o OUT --no_smk

  # 只从 accession list 获取 metadata
  python {script_name} -i accessions.tsv -o OUT --mode metadata-only

  # 已有 BLAST TSV，只做 TopN
  python {script_name} -i blast.tsv -o OUT --mode topN-only --topN 10

核心参数:
  -i, --input FILE/DIR        输入 FASTA / BLAST TSV / accession list
  -o, --output DIR            输出目录，默认 BLAST_results
  --mode MODE                 默认 blast2topN2metadata
  -bm, --blast_method METHOD  online 或 local，默认 online
  --online_engine ENGINE      默认 blast_remote
  --blast_db_nt DB            本地核酸库
  --blast_db_nr DB            本地蛋白库
  --topN N                    每个 query 保留前 N 个 hit，默认 10
  --resume                    断点续跑
  --no_smk                    不用 Snakemake
  --email EMAIL               NCBI metadata/qblast 邮箱
  --no_ncbi_metadata          不联网获取 NCBI metadata
  --local_meta_tsv FILE       使用本地 metadata TSV

说明:
  默认就是 online/remote 流程；如果要本地库，加 -bm local。
  默认启用 Snakemake 调度；如果不想用，加 --no_smk。
  查看全部高级参数:
    python {script_name} --help
"""


def parse_arguments(argv: Optional[List[str]] = None) -> argparse.Namespace:
    argv = list(sys.argv[1:] if argv is None else argv)
    script_name = Path(sys.argv[0]).name

    # Codex 修改(help-short): -h 输出短帮助，--help 保留完整 argparse 帮助。
    # 理由: 用户日常只需要核心命令；高级参数列表很长，继续放在 --help。
    # 原代码摘要: 原来 argparse 默认让 -h 和 --help 显示同一份完整帮助。
    if "-h" in argv and "--help" not in argv:
        print(short_help_text(script_name))
        raise SystemExit(0)

    # Codex 修改(help-short): 无参数运行时显示短帮助后按缺少输入报错退出。
    # 理由: 比完整 argparse usage 更易读，同时保留 CLI 参数错误的退出码 2。
    # 原代码摘要: 原来无参数直接由 argparse 报 required -i/--input 错误。
    if not argv:
        print(short_help_text(script_name), file=sys.stderr)
        print("error: the following arguments are required: -i/--input", file=sys.stderr)
        raise SystemExit(2)

    p = argparse.ArgumentParser(
        description=(
            "Batch BLAST with TopN + metadata merge, local-meta enabled.\n\n"
            "Key notes for local NT-like DB:\n"
            "  1) local_nt_metadata.tsv usually uses accession as key\n"
            "  2) sequence_metadata.tsv usually uses pseudo_accession as key\n"
            "  3) default --local_meta_key auto detects the right column\n"
            "  4) Snakemake scheduling is ON by default for safe global resource control"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "Examples\n"
            "--------\n"
            "Pure local DB:\n"
            "  python blast2_metadatav2.12.py -i query.fasta -o OUT \\\n"
            "    --mode blast2topN2metadata -bm local --blast_db_nt my_local_db \\\n"
            "    --no_ncbi_metadata --local_meta_tsv ref_out/local_nt_metadata.tsv\n\n"
            "Mixed nt + local alias:\n"
            "  python blast2_metadatav2.12.py -i query.fasta -o OUT \\\n"
            "    --mode blast2topN2metadata -bm local --blast_db_nt nt_plus_local \\\n"
            "    --local_meta_tsv ref_out/local_nt_metadata.tsv --pseudo_acc_prefixes LOC\n\n"
            "Default global scheduling (auto CPU/RAM):\n"
            "  python blast2_metadatav2.12.py -i fasta_dir -o OUT\n\n"
            "Shared-server safe mode with busy-guard:\n"
            "  python blast2_metadatav2.12.py -i fasta_dir -o OUT \\\n"
            "    --threads 64 --smk_busy_guard --smk_auto_min_mem_mb_per_job 32000\n"
        ),
    )

    p.add_argument("-i", "--input",  required=True,
                   # Codex 修改: help 补充 topN-only 可直接接收 BLAST TSV。
                   # 理由: 修复后的行为覆盖旧脚本分步使用场景，帮助信息需要同步。
                   help="Input FASTA file/directory, BLAST TSV (topN-only), or accession list (metadata-only).")
    p.add_argument("-o", "--output", default="BLAST_results",
                   help="Output directory (default: BLAST_results)")
    p.add_argument("--resume", action="store_true", help="Resume from previous run")

    p.add_argument(
        "--mode",
        choices=["blast-only", "blast2topN", "blast2meta", "blast2topN2metadata",
                 "metadata-only", "topN-only"],
        default="blast2topN2metadata",
        help="Operation mode",
    )

    p.add_argument("--extension_fasta", default=".fasta,.fas,.faa,.fa,.fna",
                   help="Comma-separated FASTA extensions (default: .fasta,.fas,.faa,.fa,.fna)")
    p.add_argument("--suffix_fasta", default="", help="Optional suffix for FASTA files")

    # BLAST
    p.add_argument("-bm", "--blast_method", choices=["local", "online"], default="online",
                   help="BLAST method (default: online)")
    p.add_argument("--online_engine", choices=["blast_remote", "biopython"],
                   default="blast_remote",
                   help="Online engine: blast_remote (BLAST+ -remote) or biopython (qblast)")
    p.add_argument("--blast_exe",    default="auto",
                   help="BLAST executable/program name (auto → blastn/blastp)")
    p.add_argument("--blast_db",     default="auto",
                   help="BLAST database name/path (auto → nt/nr)")
    p.add_argument("--blast_db_nt",  default="nt",
                   help="Local nucleotide BLAST database path/name (default: nt)")
    p.add_argument("--blast_db_nr",  default="nr",
                   help="Local protein BLAST database path/name (default: nr)")
    p.add_argument("--evalue",           default=1e-5,  type=float, help="E-value threshold")
    p.add_argument("--max_target_seqs",  default=100,   type=int,   help="Max target seqs per query")
    p.add_argument("--outfmt_blast",     default="6",
                   help='Output format: "6" or "7" uses built-in 14 columns; '
                        'custom: "6 qseqid sseqid …" must include all 14 standard fields')
    p.add_argument("-epb", "--extra_para2_blast", default="",
                   help="Extra parameters appended to BLAST command (do NOT include -outfmt)")

    # Resources
    p.add_argument("--threads", default=0, type=int,
                   help="Total CPU thread ceiling. 0 = auto-detect host CPUs")
    p.add_argument("--parallel_jobs", default=None, type=int,
                   help="Concurrent BLAST jobs inside one script run")
    p.add_argument("--local_threads_per_job", default=None, type=int,
                   help="Override per-job -num_threads for local BLAST")
    p.add_argument("--batch_size",            default=50, type=int,
                   help="Sequences per BLAST batch FASTA (default: 50)")
    p.add_argument("--online_max_concurrent", default=1,  type=int,
                   help="Max concurrent online BLAST batches (default: 1)")
    p.add_argument("--online_delay_sec",      default=12, type=int,
                   help="Delay in seconds after each online batch request (default: 12)")
    p.add_argument("--retries",    default=3,  type=int,
                   help="Retries for failed BLAST/metadata fetch (default: 3)")
    p.add_argument("--retry_sleep", default=10, type=int,
                   help="Initial retry sleep in seconds (default: 10)")

    # TopN
    p.add_argument("--topN",    default=10, type=int,
                   help="Top N hits per query (default: 10)")
    p.add_argument("--sort_by", choices=["identity", "evalue", "coverage"], default="identity",
                   help="TopN sorting criterion (default: identity)")

    # Entrez
    p.add_argument("--email",       default="you@example.com",
                   help="Email address for NCBI Entrez / qblast")
    p.add_argument("--api_key",     default="",
                   help="NCBI API key (enables 10 req/sec for Entrez E-utilities)")
    p.add_argument("--entrez_tool", default="blast2_metadatav2.12",
                   help="Entrez.tool string")
    p.add_argument("--entrez_delay_sec", default=1.0, type=float,
                   help="Extra delay between Entrez batches in seconds (default: 1.0)")

    p.add_argument("-mty", "--metadata_type", choices=["all", "short"], default="all",
                   help="Metadata verbosity: all or short (default: all)")
    p.add_argument("-sm",  "--shorter_metadata",
                   default="organism,reference_titles,reference_authors,reference_journals",
                   help="Fields included when --metadata_type short")

    # Local meta
    p.add_argument("--no_ncbi_metadata", action="store_true",
                   help="Skip Entrez efetch entirely")
    p.add_argument("--local_meta_tsv", default="",
                   help="Path to local metadata TSV (e.g. local_nt_metadata.tsv)")
    p.add_argument("--local_meta_key", default="auto",
                   help="Key column in local meta TSV. auto → accession > pseudo_accession > first col")
    p.add_argument("--local_meta_mode", choices=["map_to_ncbi", "prefix_only"],
                   default="map_to_ncbi",
                   help="How to merge local meta: map_to_ncbi or prefix_only")
    p.add_argument("--local_meta_prefix", default="local_",
                   help="Column prefix when --local_meta_mode=prefix_only")
    p.add_argument("--pseudo_acc_prefixes", default="LOC",
                   help="Comma-separated accession prefixes treated as local/pseudo")

    # Snakemake, default ON
    p.add_argument("--smk", "--snakemake", dest="use_snakemake", action="store_true",
                   help="Enable built-in Snakemake scheduling (default: ON)")
    p.add_argument("--no_smk", dest="use_snakemake", action="store_false",
                   help="Disable Snakemake and run plain Python loop")

    p.add_argument("--smk_cores", type=int, default=0,
                   help="Global Snakemake cores ceiling. 0 = auto")
    p.add_argument("--smk_jobs",  type=int, default=0,
                   help="Global concurrent Snakemake jobs. 0 = auto")
    p.add_argument("--smk_threads_per_job", type=int, default=0,
                   help="Threads per rule. 0 = auto")
    p.add_argument("--smk_total_mem_mb",  type=int, default=0,
                   help="Global mem_mb ceiling. 0 = auto")
    p.add_argument("--smk_mem_mb_per_job", type=int, default=0,
                   help="Per-job mem_mb. 0 = auto")

    p.add_argument("--smk_auto_mem", dest="smk_auto_mem", action="store_true",
                   help="Auto-detect RAM for Snakemake mem_mb budget (default: ON)")
    p.add_argument("--no_smk_auto_mem", dest="smk_auto_mem", action="store_false",
                   help="Disable automatic RAM detection")
    p.add_argument("--smk_mem_fraction", type=float, default=0.90,
                   help="Fraction of total RAM considered (default: 0.90)")
    p.add_argument("--smk_max_available_mem_fraction", type=float, default=0.85,
                   help="Max fraction of currently available RAM exposed (default: 0.85)")
    p.add_argument("--smk_auto_min_threads_per_job", type=int, default=4,
                   help="Min threads kept per job in auto mode (default: 4)")
    p.add_argument("--smk_auto_min_mem_mb_per_job", type=int, default=24576,
                   help="Min memory per job in auto mode in MB (default: 24576)")
    p.add_argument("--smk_auto_max_jobs", type=int, default=8,
                   help="Auto mode upper bound for concurrent jobs (default: 8)")

    p.add_argument("--smk_total_job_slots",  type=int, default=0,
                   help="Abstract job_slots resource total. 0 = smk_jobs")
    p.add_argument("--smk_job_slots_per_job", type=int, default=1,
                   help="job_slots consumed per job (default: 1)")
    p.add_argument("--smk_latency_wait",  type=int, default=60,
                   help="Snakemake --latency-wait in seconds (default: 60)")
    p.add_argument("--smk_restart_times", type=int, default=1,
                   help="Snakemake restart/retry count for failed rules (default: 1)")

    p.add_argument("--smk_busy_guard", dest="smk_busy_guard", action="store_true",
                   help="Auto-throttle CPU/memory when server is busy (default: ON)")
    p.add_argument("--no_smk_busy_guard", dest="smk_busy_guard", action="store_false",
                   help="Disable busy-server throttling")
    p.add_argument("--smk_busy_guard_min_cores",  type=int,   default=1,
                   help="Busy guard: minimum global cores to keep (default: 1)")
    p.add_argument("--smk_busy_guard_min_mem_mb", type=int,   default=4096,
                   help="Busy guard: minimum global memory in MB to keep (default: 4096)")
    p.add_argument("--smk_idle_cpu_fraction",       type=float, default=0.90,
                   help="Idle server: max fraction of CPU to expose (default: 0.90)")
    p.add_argument("--smk_idle_avail_mem_fraction", type=float, default=0.85,
                   help="Idle server: max fraction of available RAM to expose (default: 0.85)")

    p.add_argument("--smk_keep_going", dest="smk_keep_going", action="store_true",
                   help="Snakemake --keep-going (default: ON)")
    p.add_argument("--no_smk_keep_going", dest="smk_keep_going", action="store_false",
                   help="Disable --keep-going")
    p.add_argument("--smk_rerun_incomplete", dest="smk_rerun_incomplete", action="store_true",
                   help="Snakemake --rerun-incomplete (default: ON)")
    p.add_argument("--no_smk_rerun_incomplete", dest="smk_rerun_incomplete", action="store_false",
                   help="Disable --rerun-incomplete")
    p.add_argument("--smk_printshellcmds", dest="smk_printshellcmds", action="store_true",
                   help="Snakemake --printshellcmds (default: ON)")
    p.add_argument("--no_smk_printshellcmds", dest="smk_printshellcmds", action="store_false",
                   help="Disable --printshellcmds")
    p.add_argument("--smk_unlock_on_lock_exception",
                   dest="smk_unlock_on_lock_exception", action="store_true",
                   help="Auto-run snakemake --unlock on stale lock (default: ON)")
    p.add_argument("--no_smk_unlock_on_lock_exception",
                   dest="smk_unlock_on_lock_exception", action="store_false",
                   help="Disable automatic --unlock recovery")
    p.add_argument("--smk_dry_run",  action="store_true", help="Snakemake --dry-run")
    p.add_argument("--smk_profile",  default=None,  help="Snakemake --profile path")
    p.add_argument("--smk_extra",    default="",    help="Extra snakemake CLI arguments")

    p.set_defaults(
        use_snakemake=True,
        smk_auto_mem=True,
        smk_keep_going=True,
        smk_rerun_incomplete=True,
        smk_printshellcmds=True,
        smk_unlock_on_lock_exception=True,
        smk_busy_guard=True,
    )

    # Hidden: set by per-stem and Snakemake child invocations
    p.add_argument("--out_stem",  default="", help=argparse.SUPPRESS)

    # Fix-1: child processes load the full serialised args from this JSON
    p.add_argument("--args_json", default="", help=argparse.SUPPRESS)

    # Fix-A/G: injected by parent into child_args.json; child uses it to scale
    # the per-process NCBI rate-limit share.  Never set by user directly.
    p.add_argument("--_smk_active_jobs", default=1, type=int, help=argparse.SUPPRESS)

    return p.parse_args(argv)


# ─────────────────────────────────────────────
# main (Fix-G: split ENTREZ_LIMITER and BLAST_LIMITER)
# ─────────────────────────────────────────────
def main() -> int:
    args = parse_arguments()

    # Fix-1: if child process, load base args from JSON then apply CLI overrides
    if args.args_json:
        json_path = Path(args.args_json)
        if not json_path.exists():
            logger.error(f"--args_json file not found: {json_path}")
            return 1
        with open(json_path, "r", encoding="utf-8") as fh:
            stored: Dict[str, object] = json.load(fh)

        # Fields the child always overrides from its own CLI
        child_override_keys = {"input", "output", "out_stem", "resume", "args_json", "use_snakemake"}
        child_overrides = {k: getattr(args, k) for k in child_override_keys if hasattr(args, k)}

        # Start from stored, apply child-specific overrides
        for k, v in stored.items():
            if k not in child_override_keys:
                setattr(args, k, v)
        for k, v in child_overrides.items():
            setattr(args, k, v)

        # Child jobs never recurse into Snakemake
        args.use_snakemake = False

    args.threads = resolve_total_threads(args)
    out_dir = ensure_dir(Path(args.output).resolve())

    # Fix-G: reinitialise BOTH limiters now that we know the real args.
    # ENTREZ_LIMITER:
    #   • api_key  → NCBI allows 10 req/sec (instead of 3) for E-utilities
    #   • _smk_active_jobs → injected by the parent process into child_args.json;
    #     each child divides the global budget evenly so 10 parallel jobs never
    #     collectively exceed 3 (or 10) req/sec.
    # BLAST_LIMITER:
    #   • STRICT 1 req / 3 sec limit, no API key support
    #   • _smk_active_jobs → divides this limit across parallel jobs
    global ENTREZ_LIMITER, BLAST_LIMITER
    active_jobs    = int(getattr(args, "_smk_active_jobs", 1))
    
    ENTREZ_LIMITER = _make_entrez_limiter(api_key=getattr(args, "api_key", ""),
                                          active_jobs=active_jobs)
    BLAST_LIMITER  = _make_blast_limiter(active_jobs=active_jobs)
    
    _entrez_max_reqs = 10 if getattr(args, "api_key", "") else 3

    fh = logging.FileHandler(out_dir / "batch_blast.log", encoding="utf-8")
    fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(fh)

    exts  = [x.strip() for x in args.extension_fasta.split(",") if x.strip()]
    files = find_input_files(args.input, exts, args.suffix_fasta, args.mode)
    if not files:
        logger.error("No input files found.")
        return 1

    preflight_checks(args, files, out_dir)

    logger.info(f"[FILES]          {len(files)}")
    logger.info(f"[PANDAS]         {HAVE_PANDAS}")
    logger.info(f"[PSUTIL]         {HAVE_PSUTIL}")
    logger.info(f"[THREAD_CEILING] {args.threads}")
    logger.info(f"[SCHEDULER]      {'snakemake' if args.use_snakemake else 'python'}")
    logger.info(
        f"[ENTREZ_RATE_LIMIT] {_entrez_max_reqs} calls / "
        f"{1.05 * active_jobs:.2f}s period "
        f"(api_key={'yes' if getattr(args,'api_key','') else 'no'}, "
        f"active_jobs={active_jobs})"
    )
    logger.info(
        f"[BLAST_RATE_LIMIT] 1 call / "
        f"{3.1 * active_jobs:.2f}s period "
        f"(strict NCBI BLAST limit, no API key support, "
        f"active_jobs={active_jobs})"
    )

    if args.use_snakemake:
        return int(run_with_snakemake(args, files, out_dir))

    dirs = setup_output_directory(out_dir)
    # Codex 修改: 记录普通 Python 模式下的文件级失败数。
    # 理由: 原代码只记录异常并继续，最后仍返回 0，自动化流程会把失败任务误判为成功。
    # 原代码:
    #   for fp in files:
    #       try:
    #           process_one_fasta(fp, args, dirs)
    #       except ...
    failed = 0
    for fp in files:
        try:
            if not process_one_fasta(fp, args, dirs):
                failed += 1
        except KeyboardInterrupt:
            logger.error("Interrupted by user.")
            return 130
        except Exception as e:
            failed += 1
            logger.error(f"[FILE_FAIL] {fp}: {e}", exc_info=True)

    if failed:
        logger.error(f"Processing completed with {failed} failed file(s).")
        return 1

    logger.info("All processing completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
