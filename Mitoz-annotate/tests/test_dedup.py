"""Tests for batch_mitoz.py deduplication (RefSeq NC_ kept over the INSDC submission it was made from).

Run: python3 -m pytest -q tests/test_dedup.py
"""
import copy
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import batch_mitoz as bm  # noqa: E402

DATA = HERE / "data"


def insdc_copy(refseq, acc: str, comment: str = ""):
    """The INSDC submission a RefSeq record was derived from (same sequence, own accession)."""
    r = copy.deepcopy(refseq)
    r.id, r.name = f"{acc}.1", acc
    r.annotations["accessions"] = [acc]
    r.annotations["comment"] = comment
    return r


def test_refseq_source_parsed():
    rec = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
    assert bm.is_refseq_acc("NC_000933") and not bm.is_refseq_acc("AF000023")
    assert bm.refseq_source_accessions(rec) == ["AF000023"]      # "derived from AF000023"


def test_refseq_kept_source_dropped_both_orders(tmp_path):
    nc = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
    src = insdc_copy(nc, "AF000023")
    other = insdc_copy(nc, "ZZ999999")
    other.seq = Seq(str(nc.seq)[::-1])                            # unrelated, different sequence
    for order in ([src, nc, other], [nc, other, src]):
        f = tmp_path / "in.gb"
        SeqIO.write(order, str(f), "genbank")
        kept, _, stats = bm.dedup_gb_records([f], dedup_report_path=tmp_path / "rep.tsv")
        ids = sorted(m.acc_versioned for m in kept)
        assert ids == ["NC_000933.1", "ZZ999999.1"], ids
        assert stats["dup_refseq"] == 1
        assert "AF000023.1 superseded by RefSeq NC_000933.1" in (tmp_path / "rep.tsv").read_text()


def test_identical_sequence_prefers_refseq_without_comment(tmp_path):
    nc = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
    nc.annotations["comment"] = ""                                 # no "derived from" link
    src = insdc_copy(nc, "AB123456")
    f = tmp_path / "in.gb"
    SeqIO.write([src, nc], str(f), "genbank")
    kept, _, stats = bm.dedup_gb_records([f])
    assert [m.acc_versioned for m in kept] == ["NC_000933.1"]
    assert stats["dup_seq"] == 1


def test_dedup_subcommand_writes_gb_and_fasta(tmp_path):
    nc = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
    f = tmp_path / "in.gb"
    SeqIO.write([insdc_copy(nc, "AF000023"), nc], str(f), "genbank")
    subprocess.run([sys.executable, "-W", "ignore", str(HERE.parent / "batch_mitoz.py"), "dedup",
                    "-i", str(f), "-o", str(tmp_path / "out")], check=True, capture_output=True)
    assert [r.id for r in SeqIO.parse(str(tmp_path / "out" / "dedup.gb"), "genbank")] == ["NC_000933.1"]
    assert [r.id for r in SeqIO.parse(str(tmp_path / "out" / "dedup.fasta"), "fasta")] == ["NC_000933.1"]
