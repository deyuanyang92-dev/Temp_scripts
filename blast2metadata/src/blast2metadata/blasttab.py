"""BLAST tabular parsing, merging and TopN ranking."""
from __future__ import annotations

import argparse
import bisect
import csv
import math
import os
import re
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Union


from .util import HAVE_PANDAS, pd
from .accession import accession_base, parse_accession_from_sseqid
from .fields import BLAST_FIELDS, BLAST_HEADERS, NUMERIC_FLOAT_COLS, NUMERIC_INT_COLS, OUTPUT_FIELD_DESC, _BLAST_LABEL_TO_CODE
from .util import ensure_dir

_HEADER_CODES = set(BLAST_FIELDS) | set(OUTPUT_FIELD_DESC) | {"metadata_length"}

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
    # 一次完成转换与校验；避免每行重复构造列集合、查找列下标及二次扫描数字列。
    for k, v in zip(col_names, cols):
        if k in NUMERIC_FLOAT_COLS or k in NUMERIC_INT_COLS:
            if v == "":
                r[k] = None
                continue
            try:
                n = float(v)
                if not math.isfinite(n) or (k in NUMERIC_INT_COLS and not n.is_integer()):
                    raise ValueError
                r[k] = int(n) if k in NUMERIC_INT_COLS else n
            except (ValueError, OverflowError):
                raise ValueError(f"invalid numeric column {k}: {v!r}") from None
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
    if len(toks) >= 2 and all(t in _HEADER_CODES for t in toks):
        return toks
    return None


def iter_blast_results(
    blast_output: Union[str, Path], col_names: Optional[List[str]] = None, quoted: bool = False,
) -> Iterator[Dict[str, object]]:
    """流式读取 BLAST TSV/CSV；坏行报文件名和行号，避免把损坏结果当作无命中。

    quoted=True 仅用于读回本软件 csv.writer 写出的表；BLAST 原始输出不转义引号，
    NCBI title 常以 "Candidatus ..." 开头，按引号解析会报错或吞行。
    """
    columns = list(col_names or BLAST_HEADERS)
    p = Path(blast_output)
    with open(p, "r", encoding="utf-8-sig", newline="") as f:
        # 修改 (2026-10-22): csv.reader 支持本软件写出的带引号、制表符和换行的 title。
        # 旧 split('\t') 不能读回 csv.writer 生成的合法 TSV。注释中的逗号不用于猜分隔符。
        delimiter = "\t"
        for line in f:
            if line.strip() and not line.startswith("#"):
                delimiter = "\t" if "\t" in line else ","
                break
        f.seek(0)
        reader = csv.reader(f, delimiter=delimiter, strict=True,
                            quoting=csv.QUOTE_MINIMAL if quoted else csv.QUOTE_NONE)
        have_schema = False
        try:
            for parts in reader:
                if not parts or not any(parts):
                    continue
                if parts[0].startswith("#"):
                    if parts[0].startswith("# Fields:"):
                        found = header_columns(delimiter.join(parts))
                        if not found:
                            raise ValueError(f"{p}: line {reader.line_num}: unknown # Fields header")
                        columns, have_schema = found, True
                    continue
                found = header_columns("\t".join(parts)) if parts[0] in _HEADER_CODES else None
                if found:
                    columns, have_schema = found, True
                    continue
                if not have_schema:
                    if len(parts) == 12 and columns == BLAST_HEADERS:
                        # 标准 outfmt 6 的 12 列没有 qcovs，不能伪造覆盖度为 0/100。
                        columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split()
                    have_schema = True
                try:
                    row = _parse_outfmt6_line(parts, columns)
                except ValueError as exc:
                    raise ValueError(f"{p}: line {reader.line_num}: {exc}") from exc
                if row is None or (len(parts) > len(columns) and "stitle" not in columns):
                    raise ValueError(f"{p}: line {reader.line_num}: expected {len(columns)} columns, got {len(parts)}")
                for target, aliases in (("qseqid", ("qaccver", "qacc")), ("sseqid", ("saccver", "sacc"))):
                    if target not in row:
                        for alias in aliases:
                            if alias in row:
                                row[target] = row[alias]
                                break
                    if not row.get(target):
                        raise ValueError(f"{p}: line {reader.line_num}: missing {target}")
                yield row
        except csv.Error as exc:
            raise ValueError(f"{p}: line {reader.line_num}: {exc}") from exc


def parse_blast_results_py(
    blast_output: Union[str, Path], col_names: Optional[List[str]] = None, quoted: bool = False,
) -> List[Dict[str, object]]:
    """兼容旧 API；TopN 管道直接消费迭代器，只有需要全表时才转 list。"""
    return list(iter_blast_results(blast_output, col_names, quoted))


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
            if p.stat().st_size == 0:
                continue
            with open(p, "r", encoding="utf-8") as f:
                for line in f:
                    if not line.strip() or line.startswith("#"):
                        continue
                    out.write(line if line.endswith("\n") else line + "\n")
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
# pident_x_qcovs 是排序启发式：代表 HSP 的 identity × subject 全部 HSP 的覆盖度。
# 多 HSP identity 不同、重叠、gap 及 qcovs 取整都会影响它，不能当作全 query 精确一致率。
DERIVED_SORT_COLS = {"pident_x_qcovs"}


def parse_sort_spec(spec: str) -> List[str]:
    """'qcovs,pident' / preset name -> list of sort columns. Raises ValueError on unknown columns."""
    spec = (spec or "").strip()
    spec = SORT_PRESETS.get(spec, spec)
    cols = [SORT_ALIASES.get(c.strip(), c.strip()) for c in spec.split(",") if c.strip()]
    if not cols:
        raise ValueError("empty --sort_by")
    known = NUMERIC_FLOAT_COLS | NUMERIC_INT_COLS | DERIVED_SORT_COLS
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
    rows: Iterable[Dict[str, object]], topN: int, sort_by: str,
    min_pident: float = 0.0, min_qcovs: float = 0.0,
) -> List[Dict[str, object]]:
    """每条 query 只保留 N 个 subject（每个 subject 取最佳 HSP 作代表，再按阈值过滤）；不缓存全部 HSP。"""
    cols = parse_sort_spec(sort_by)
    if topN <= 0:
        raise ValueError("topN must be positive")
    required = set(cols) - DERIVED_SORT_COLS
    if "pident_x_qcovs" in cols:
        required.update(("pident", "qcovs"))
    if min_pident > 0:
        required.add("pident")
    if min_qcovs > 0:
        required.add("qcovs")

    def number(v):
        try:
            n = float(v)
            return n if math.isfinite(n) else None
        except (ValueError, TypeError, OverflowError):
            return None

    def ranked():
        for original in rows:
            missing = required - original.keys()
            if missing:
                raise ValueError(f"TopN requires columns {sorted(missing)}; supply them in the BLAST table")
            row = original
            if "pident_x_qcovs" in cols:
                pid, cov = number(original.get("pident")), number(original.get("qcovs"))
                row = dict(original)
                row["pident_x_qcovs"] = round(pid * cov / 100, 3) if pid is not None and cov is not None else None
            rank = []
            for col in cols:
                v = number(row.get(col))
                rank.extend((v is None, 0.0 if v is None else v if col in SORT_ASCENDING else -v))
            q = str(row.get("qseqid", ""))
            subject = str(row.get("sseqid", ""))
            if not q or not subject:
                raise ValueError("TopN requires non-empty qseqid and sseqid")
            # 精确 subject ID 去重：去掉 . 后缀会误合并 accession 版本及 contig.partA/partB。
            # 同分时按 ID/坐标/title 固定次序，避免并行批次完成顺序影响最终命中。
            key = tuple(rank) + (subject,) + tuple(str(row.get(c, "")) for c in
                                                  ("qstart", "qend", "sstart", "send", "stitle"))
            yield q, subject, key, row

    # 修改 (2026-10-23): 先选每个 query-subject 的代表 HSP，再用 --min_pident/--min_qcovs 过滤。
    # 理由: 原先逐 HSP 过滤后再选代表。subject 的主 HSP 900 bp/85% 未达 97% 时，旁边一段
    #       40 bp/100% 的短 HSP 仍能通过；qcovs 是 subject 级（含 85% 那段）也能通过，
    #       于是在“库里没有同种序列”时报出假的高一致率命中。现与 NCBI 网页一致：
    #       subject 的 Per. Ident 取其最佳 HSP。
    # BLAST+/DIAMOND/MMseqs2 的表格把同一 query-subject 的 HSP 连续输出，所以只需合并
    # 相邻行，内存仍为 O(query数 × N)。行顺序被打乱的外部表会分成多段，每段各自选代表，
    # 此时退化为原先的逐段过滤，不会比原来更差。
    def representatives():
        current = None
        for item in ranked():
            if current is not None and item[:2] == current[:2]:
                if item[2] < current[2]:
                    current = item
                continue
            if current is not None:
                yield current
            current = item
        if current is not None:
            yield current

    # 修改 (2026-10-22): 每个 query 维护有序的 N 个候选和 subject 索引。
    # 已淘汰 subject 的旧 HSP 不可能超过随后只会提高的门槛；它的新 HSP 仍正常竞争。
    # 内存 O(query数 × N)，不再为所有 HSP 建字典列表并全量排序；插入 O(N)，适合小 TopN。
    groups = {}
    for q, subject, key, row in representatives():
        pid, cov = number(row.get("pident")), number(row.get("qcovs"))
        if min_pident > 0 and (pid is None or pid < min_pident):
            continue
        if min_qcovs > 0 and (cov is None or cov < min_qcovs):
            continue
        ordered, selected = groups.setdefault(q, ([], {}))
        previous = selected.get(subject)
        if previous is not None:
            if key >= previous[0]:
                continue
            ordered.remove((previous[0], subject))
        elif len(ordered) == topN:
            if key >= ordered[-1][0]:
                continue
            _, discarded = ordered.pop()
            del selected[discarded]
        bisect.insort(ordered, (key, subject))
        selected[subject] = (key, row)
    return [selected[s][1] for q in sorted(groups)
            for ordered, selected in [groups[q]] for _, s in ordered]


def extract_accession_numbers_from_rows(rows: List[Dict[str, object]]) -> List[str]:
    clean: List[str] = []
    seen = set()  # O(1) membership instead of quadratic list scans on large hit tables.
    for r in rows:
        acc  = parse_accession_from_sseqid(str(r.get("sseqid", "")))
        a    = accession_base(acc)
        if a and a not in seen:
            seen.add(a)
            clean.append(a)
    return clean
