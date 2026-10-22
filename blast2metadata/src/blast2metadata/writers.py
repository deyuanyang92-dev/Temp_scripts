"""Delimited output writers."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Optional, Sequence


from .accession import accession_base, is_local_or_pseudo_accession, parse_accession_from_sseqid
from .fields import BLAST_HEADERS, NT_META_COLS, column_label
from .util import ensure_dir, logger

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
    header_style: str = "code",
    encoding: str = "utf-8",
) -> None:
    ensure_dir(path.parent)
    ordered = _ordered_row_keys(rows, default_headers=default_headers)
    with open(path, "w", encoding=encoding, newline="") as f:
        w = csv.writer(f, delimiter=delimiter)
        meta_table = list(default_headers or []) == NT_META_COLS
        w.writerow([column_label(k, header_style, meta_table) for k in ordered])
        for r in rows:
            w.writerow(["" if r.get(k) is None else r.get(k) for k in ordered])


def write_results_csv(
    path: Path,
    rows: List[Dict[str, object]],
    args: argparse.Namespace,
    default_headers: Optional[Sequence[str]] = None,
) -> None:
    """Final *_results.csv. 修改 (2026-10-11): 统一写出并按 --header_style 生成表头；
    utf-8-sig 让 Excel 直接双击打开时正确识别编码。原代码 5 处分别用 pandas.to_csv / rows_to_csv。"""
    rows_to_delimited(path, rows, ",", default_headers=default_headers,
                      header_style=getattr(args, "header_style", "code"), encoding="utf-8-sig")


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
