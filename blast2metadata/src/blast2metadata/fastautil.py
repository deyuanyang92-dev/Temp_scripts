"""FASTA helpers: sequence-type detection and batching."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List, Union

from Bio import SeqIO  # type: ignore

from .util import ensure_dir, logger, open_text

def detect_sequence_type(fasta_file: Union[str, Path], max_records: int = 2) -> str:
    fasta_file = str(fasta_file)
    try:
        with open_text(fasta_file) as f:
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
        base    = fasta_file.name[:-3] if fasta_file.name.lower().endswith(".gz") else fasta_file.name
        out_fa  = temp_dir / f"{bid}__{base}"
        with open(out_fa, "w", encoding="utf-8") as fh:
            SeqIO.write(cur_records, fh, "fasta")
        batches.append(BatchInfo(batch_id=bid, fasta_path=out_fa))
        bidx += 1
        cur_records = []

    # 修改 (2026-10-22): 支持 .gz 输入；重复 query ID 会让 TopN 把不同序列的命中混在一起，提前警告。
    seen_ids: set = set()
    dup_ids: List[str] = []
    with open_text(fasta_file) as f:
        for rec in SeqIO.parse(f, "fasta"):
            if rec.id in seen_ids:
                dup_ids.append(rec.id)
            seen_ids.add(rec.id)
            cur_records.append(rec)
            if len(cur_records) >= batch_size:
                flush()
        flush()
    if dup_ids:
        logger.warning(
            f"[FASTA] {fasta_file.name}: {len(dup_ids)} duplicate sequence ID(s) "
            f"(e.g. {dup_ids[:3]}); hits of same-ID sequences are pooled in TopN"
        )

    return batches
