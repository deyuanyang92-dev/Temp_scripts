"""Orientation: fallback genes when cox1 is absent, and the 5'-most piece of a split gene.

Run: python3 -m pytest -q tests/test_orientation.py
"""
import argparse
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import batch_mitoz as bm  # noqa: E402


def rec_with(features):
    r = SeqRecord(Seq("ACGT" * 500), id="R1", name="R1", description="test")
    r.annotations.update(molecule_type="DNA", topology="circular")
    r.features = [SeqFeature(FeatureLocation(a, b, strand=st), type=t, qualifiers={"gene": [g]})
                  for t, g, a, b, st in features]
    return r


def test_aliases():
    assert "coii" in bm.gene_aliases("cox2") and "cytb" in bm.gene_aliases("cob")
    assert bm.gene_aliases("COX1") == bm.gene_aliases("COI")


def test_orientation_genes_from_args():
    a = argparse.Namespace(which_gene_first="cox1", fallback_genes="cox2,cob")
    assert bm.orientation_genes(a) == "cox1,cox2,cob"
    assert bm.orientation_genes(argparse.Namespace(which_gene_first="cox1", fallback_genes="none")) == "cox1"


def test_fallback_when_cox1_absent():
    r = rec_with([("CDS", "COII", 300, 900, 1), ("CDS", "CYTB", 1000, 1500, 1)])
    feat, gene = bm.find_orientation_feature(r, "cox1,cox2,cob")
    assert gene == "cox2" and int(feat.location.start) == 300


def test_split_gene_takes_five_prime_piece():
    plus = rec_with([("CDS", "cox1", 1200, 1700, 1), ("CDS", "cox1", 200, 700, 1)])   # exon 2 listed first
    assert int(bm.find_orientation_feature(plus, "cox1")[0].location.start) == 200
    minus = rec_with([("CDS", "cox1", 200, 700, -1), ("CDS", "cox1", 1200, 1700, -1)])
    assert int(bm.find_orientation_feature(minus, "cox1")[0].location.end) == 1700


def test_nothing_found():
    assert bm.find_orientation_feature(rec_with([("CDS", "nad5", 10, 400, 1)]), "cox1,cox2") == (None, "")


def test_reorient_uses_fallback(tmp_path):
    bm.ORIENTATION_FALLBACK.clear()
    r = rec_with([("CDS", "cox2", 300, 900, 1), ("CDS", "nad5", 1000, 1500, 1)])
    gb_in = tmp_path / "in.gbf"
    SeqIO.write(r, str(gb_in), "genbank")
    n_ok, n_no, _ = bm.reorient_one_gb(gb_in, tmp_path / "out.gbf", tmp_path / "no.gbf", "cox1,cox2", True)
    assert (n_ok, n_no) == (1, 0)
    out = SeqIO.read(str(tmp_path / "out.gbf"), "genbank")
    assert str(out.seq)[:20] == str(r.seq)[300:320]               # sequence now starts at cox2
    assert bm.ORIENTATION_FALLBACK == {"R1": "cox2"}
