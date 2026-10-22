"""`blast2metadata --doctor`: report external tools and Python packages found on this machine."""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from typing import List, Optional, Tuple

from . import __version__

# (name, [command to print version], regex to extract version)
_TOOLS: List[Tuple[str, List[str], str]] = [
    ("blastn",      ["blastn", "-version"],      r"blastn:\s*(\S+)"),
    ("blastp",      ["blastp", "-version"],      r"blastp:\s*(\S+)"),
    ("makeblastdb", ["makeblastdb", "-version"], r"makeblastdb:\s*(\S+)"),
    ("diamond",     ["diamond", "version"],      r"version\s+(\S+)"),
    ("mmseqs",      ["mmseqs", "version"],       r"(\S+)"),
    ("snakemake",   ["snakemake", "--version"],  r"(\S+)"),
]


def _tool_version(cmd: List[str], pat: str) -> Optional[str]:
    exe = shutil.which(cmd[0])
    if not exe:
        return None
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        m = re.search(pat, (p.stdout or "") + (p.stderr or ""))
        return (m.group(1) if m else "found") + f"  ({exe})"
    except Exception as e:  # pragma: no cover - defensive
        return f"found but `{' '.join(cmd)}` failed: {e}  ({exe})"


def _py_version(mod: str) -> Optional[str]:
    try:
        m = __import__(mod)
    except Exception:
        return None
    return str(getattr(m, "__version__", "found"))


def run_doctor(exes: Optional[dict] = None) -> int:
    """exes: {'diamond': path, 'mmseqs': path} overrides from --diamond_exe / --mmseqs_exe."""
    exes = exes or {}
    print(f"blast2metadata {__version__}  (python {sys.version.split()[0]}, {sys.executable})")
    print("\nExternal tools:")
    for name, cmd, pat in _TOOLS:
        cmd = [exes.get(name, cmd[0])] + cmd[1:]
        v = _tool_version(cmd, pat)
        note = ""
        if name == "diamond" and v:
            m = re.match(r"(\d+)\.(\d+)\.(\d+)", v)
            if m:
                ok = tuple(int(x) for x in m.groups()) >= (2, 1, 14)
                note = ("  [reads BLAST v5 protein DBs directly: --blast_db /path/nr]" if ok
                        else "  [< 2.1.14: cannot read BLAST DBs, needs .dmnd]")
        print(f"  {name:12s} {v if v else 'NOT FOUND'}{note}")
    print("\nPython packages:")
    for mod, note in [("Bio", "biopython, required"), ("pandas", "optional"),
                      ("psutil", "optional, host load detection"), ("snakemake", "optional scheduler")]:
        v = _py_version(mod)
        print(f"  {mod:12s} {v if v else 'NOT FOUND'}  [{note}]")
    print("\nEnvironment:")
    proxies = {k: os.environ[k] for k in ("https_proxy", "HTTPS_PROXY", "http_proxy", "HTTP_PROXY") if os.environ.get(k)}
    print(f"  HTTP proxy   {proxies if proxies else 'none'}"
          + ("  (blastn -remote may fail behind a proxy; use --online_engine biopython)" if proxies else ""))
    print(f"  BLASTDB      {os.environ.get('BLASTDB', 'not set')}")
    print(f"  NCBI_API_KEY {'set' if os.environ.get('NCBI_API_KEY') else 'not set'}")
    return 0
