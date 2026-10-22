"""Accession parsing, accession-list input and Entrez DB guessing."""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union


from .fields import FASTA_EXTENSIONS_FALLBACK, NUCLEOTIDE_ACCESSION_PREFIXES, PROTEIN_ACCESSION_PREFIXES

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
    name = str(path)
    if name.lower().endswith(".gz"):
        name = name[:-3]
    return path_has_extension(name, extensions or FASTA_EXTENSIONS_FALLBACK)


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
