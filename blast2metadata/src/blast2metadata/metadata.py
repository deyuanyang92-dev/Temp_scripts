"""Local metadata TSV, Entrez GenBank metadata fetch and merge."""
from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from Bio import Entrez  # type: ignore

from . import ratelimit
from .accession import accession_base, guess_entrez_db_for_accession, parse_accession_from_sseqid
from .engines import entrez_setup
from .fields import NT_META_COLS
from .util import logger

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
    batch_size: int = 100,
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
        ratelimit.ENTREZ_LIMITER.acquire()
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
        # 修改 (2026-10-22): 整批失败（常见: 一个无效 accession 让 efetch 返回 HTTP 400）时对半拆分重试，
        # 原代码直接放弃整批，其余有效 accession 也拿不到 metadata。
        if len(ids) > 1:
            half = len(ids) // 2
            logger.info(f"[META_FETCH] {label} db={db_name}: splitting {len(ids)} ids and retrying")
            left, _ = fetch_ids_from_db(ids[:half], db_name, label + ".a")
            right, _ = fetch_ids_from_db(ids[half:], db_name, label + ".b")
            if left or right:
                return left + right, None
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
