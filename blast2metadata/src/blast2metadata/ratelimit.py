"""NCBI rate limiters (Entrez E-utilities vs. BLAST web API)."""
from __future__ import annotations

import time
import threading
from typing import List


# ─────────────────────────────────────────────
# Fix-G: Split NCBI rate limiter into ENTREZ and BLAST
# ─────────────────────────────────────────────
# The NCBI API Key system *only applies to Entrez E-utilities* (esearch, efetch).
# The NCBI BLAST API (https://blast.ncbi.nlm.nih.gov/Blast.cgi) is a completely
# separate backend that does NOT support API keys and enforces a strict limit
# of 1 request every 10 seconds. Sharing a single rate limiter between them
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
    
    CRITICAL: The BLAST URL API allows at most 1 request every 10 seconds,
    regardless of API keys. To be safe across distributed clock jitter, we use 
    10.1 seconds per request.
    
    Parameters
    ----------
    active_jobs : number of concurrently running Snakemake child processes;
                  each gets 1 request per (10.1 × active_jobs) seconds.
    """
    active_jobs = max(1, int(active_jobs))
    # 修改 (2026-10-23): 3.1 s → 10.1 s。
    # 理由: NCBI BLAST URL API 开发者说明要求 "Do not contact the server more often than
    #       once every 10 seconds"（Biopython qblast 源码注释同样引用此条）。原来 3.1 s 只因
    #       每次请求后另有 --online_delay_sec 12 s 才没有超速；现在去掉了那段重复的固定等待，
    #       由这里单独保证提交间隔。轮询由 Biopython / BLAST+ 自己节流。
    safe_period = 10.1 * active_jobs
    return _RateLimiter(max_calls=1, period=safe_period)


# Module-level sentinels — replaced in main() after args are parsed.
ENTREZ_LIMITER: _RateLimiter = _make_entrez_limiter()
BLAST_LIMITER: _RateLimiter = _make_blast_limiter()
