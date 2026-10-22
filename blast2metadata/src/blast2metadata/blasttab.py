"""BLAST tabular parsing, merging and TopN ranking."""
from __future__ import annotations

import argparse
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union


from .util import HAVE_PANDAS, pd
from .accession import accession_base, parse_accession_from_sseqid
from .fields import BLAST_FIELDS, BLAST_HEADERS, NUMERIC_FLOAT_COLS, NUMERIC_INT_COLS, OUTPUT_FIELD_DESC, _BLAST_LABEL_TO_CODE
from .util import _to_float, _to_int, ensure_dir

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


def header_columns(line: str) -> Optional[List[str]]:
    """Column codes if *line* is a header: outfmt 7 '# Fields:' or a row of field codes."""
    s = line.strip().lstrip("\ufeff")
    if s.startswith("# Fields:"):
        labels = [x.strip() for x in s[len("# Fields:"):].split(",")]
        codes = [_BLAST_LABEL_TO_CODE.get(x) for x in labels]
        return None if None in codes else [c for c in codes if c]
    if not s or s.startswith("#"):
        return None
    toks = [t.strip() for t in re.split(r"[\t,]", s)]
    known = set(BLAST_FIELDS) | set(OUTPUT_FIELD_DESC) | {"metadata_length"}
    if len(toks) >= 2 and all(t in known for t in toks):
        return toks
    return None


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
    with open(p, "r", encoding="utf-8-sig", errors="ignore") as f:
        for line in f:
            # 修改 (2026-10-11): 识别表头行和 outfmt 7 的 "# Fields:" 行，按其列名解析。
            # 理由: 合并表现在带表头；--resume 读回 TopN TSV 时原代码把表头当成一条数据。
            fields = header_columns(line)
            if fields:
                if not rows:
                    col_names = fields
                continue
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


def merge_blast_tsv_text(
    blast_files: Sequence[Union[str, Path]], merged_out: Path, header: Optional[Sequence[str]] = None
) -> None:
    ensure_dir(merged_out.parent)
    part = merged_out.with_name(merged_out.name + ".part")
    with open(part, "w", encoding="utf-8") as out:
        if header:
            out.write("\t".join(header) + "\n")
        for fp in blast_files:
            p = Path(fp)
            if not p.exists() or p.stat().st_size == 0:
                continue
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if not line.strip() or line.startswith("#"):
                        continue
                    out.write(line)
    os.replace(part, merged_out)


# ─────────────────────────────────────────────
# TopN
# ─────────────────────────────────────────────
# 修改 (2026-10-11): TopN 排序可配置。--sort_by 接受预设名或逗号分隔的列（按顺序作主次排序键），
# 方向自动：evalue/mismatch/gapopen/gaps 越小越好，其余越大越好。
# 预设与原 4 种排序逐键一致（identity/evalue/coverage/bitscore）。
SORT_PRESETS: Dict[str, str] = {
    "bitscore": "bitscore,evalue,pident",
    "identity": "pident,qcovs,evalue",
    "coverage": "qcovs,pident,evalue",
    "evalue":   "evalue,pident,qcovs",
    "pident_x_qcovs": "pident_x_qcovs,bitscore,evalue",
}
SORT_ALIASES: Dict[str, str] = {"identity": "pident", "coverage": "qcovs"}
SORT_ASCENDING = {"evalue", "mismatch", "gapopen", "gaps"}
# pident_x_qcovs = pident * qcovs / 100: % of the query that is identical to the subject
# (identity and coverage in one number). Exact for single-HSP hits; for multi-HSP hits it
# uses the best HSP's pident with the per-subject coverage.
DERIVED_SORT_COLS = {"pident_x_qcovs"}


def parse_sort_spec(spec: str) -> List[str]:
    """'qcovs,pident' / preset name -> list of sort columns. Raises ValueError on unknown columns."""
    spec = (spec or "").strip()
    spec = SORT_PRESETS.get(spec, spec)
    cols = [SORT_ALIASES.get(c.strip(), c.strip()) for c in spec.split(",") if c.strip()]
    if not cols:
        raise ValueError("empty --sort_by")
    known = set(BLAST_FIELDS) | DERIVED_SORT_COLS
    bad = [c for c in cols if c not in known]
    if bad:
        raise ValueError(
            f"unknown --sort_by column(s) {bad}; use a preset ({', '.join(SORT_PRESETS)}) "
            "or BLAST columns, e.g. qcovs,pident"
        )
    return cols


def _sort_spec_arg(spec: str) -> str:
    try:
        parse_sort_spec(spec)
    except ValueError as e:
        raise argparse.ArgumentTypeError(str(e))
    return spec


def get_topN_hits_py(
    rows: List[Dict[str, object]], topN: int, sort_by: str,
    min_pident: float = 0.0, min_qcovs: float = 0.0,
) -> List[Dict[str, object]]:
    if not rows:
        return rows
    cols = parse_sort_spec(sort_by)

    def fnum(v, default: float) -> float:
        if v is None:
            return default
        try:
            return float(v)
        except Exception:
            return default

    if min_pident > 0 or min_qcovs > 0:
        rows = [r for r in rows
                if fnum(r.get("pident"), -1) >= min_pident and fnum(r.get("qcovs"), -1) >= min_qcovs]
    if "pident_x_qcovs" in cols:
        rows = [dict(r) for r in rows]
        for r in rows:
            pid, cov = fnum(r.get("pident"), -1), fnum(r.get("qcovs"), -1)
            r["pident_x_qcovs"] = round(pid * cov / 100, 3) if pid >= 0 and cov >= 0 else None

    # missing values sort last whatever the direction
    def sort_key(r):
        k: List[object] = [str(r.get("qseqid", ""))]
        for c in cols:
            v = fnum(r.get(c), float("nan"))
            missing = v != v
            k += [missing, 0.0 if missing else (v if c in SORT_ASCENDING else -v)]
        return tuple(k)

    sorted_rows = sorted(rows, key=sort_key)

    # 修改 (2026-10-10): TopN 按 subject 计数，每个 query-subject 只保留排序最优的一个 HSP。
    # 理由: outfmt 6 每行是一个 HSP；原代码按行计数，跨环形起点或含内含子的命中会拆成多个 HSP，
    #       top N 被同一 subject 占满（实测 LC670763 COX1 top3 中同一 subject 出现两次）。
    # 原代码:
    #   for r in sorted_rows:
    #       q = str(r.get("qseqid", "")); c = seen.get(q, 0)
    #       if c < topN: out.append(r); seen[q] = c + 1
    out: List[Dict[str, object]] = []
    seen: Dict[str, int] = {}
    seen_pairs: set = set()
    for r in sorted_rows:
        q = str(r.get("qseqid", ""))
        pair = (q, accession_base(parse_accession_from_sseqid(str(r.get("sseqid", "")))))
        if pair in seen_pairs:
            continue
        c = seen.get(q, 0)
        if c < topN:
            out.append(r)
            seen[q] = c + 1
            seen_pairs.add(pair)
    return out


def extract_accession_numbers_from_rows(rows: List[Dict[str, object]]) -> List[str]:
    clean: List[str] = []
    for r in rows:
        acc  = parse_accession_from_sseqid(str(r.get("sseqid", "")))
        a    = accession_base(acc)
        if a and a not in clean:
            clean.append(a)
    return clean
