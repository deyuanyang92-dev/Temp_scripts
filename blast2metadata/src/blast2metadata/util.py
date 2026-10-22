"""Logging, subprocess helper, host-state detection and small helpers."""
from __future__ import annotations

import argparse
import gzip
import logging
import os
import shutil
import subprocess
import sys
import time
import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union


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

# ─────────────────────────────────────────────
# logging
# ─────────────────────────────────────────────
logger = logging.getLogger("blast2_metadata")
logger.setLevel(logging.INFO)
_handler = logging.StreamHandler(sys.stdout)
_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logger.handlers.clear()
logger.addHandler(_handler)
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


def open_text(path: Union[str, Path], mode: str = "r"):
    """Open a text file, transparently decompressing *.gz (FASTA inputs may be gzipped)."""
    p = str(path)
    if p.lower().endswith(".gz"):
        return gzip.open(p, mode + "t", encoding="utf-8", errors="ignore")
    return open(p, mode, encoding="utf-8", errors="ignore")


def input_stem(path: Union[str, Path]) -> str:
    """'x.fasta.gz' -> 'x'; 'x.fasta' -> 'x' (output names never carry the .gz)."""
    name = Path(path).name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    return Path(name).stem


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
        # 修改 (2026-10-10): 每次失败都记 stderr 末尾。理由: 原代码只记 rc，blastn -remote rc=255
        # 时看不到原因（实测为代理下 "Connection stream is in bad state"）。
        err_tail = (p.stderr or "").strip()[-500:]
        logger.warning(f"[CMD_FAIL] attempt {attempt}/{retries} rc={p.returncode} {err_tail}")
        if attempt < retries:
            time.sleep(sleep)
            sleep = int(sleep * 1.8)
    return last


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
def _write_text(path: Path, text: str) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
