"""Tests for mito_rescue.py (genes MitoZ missed, added from MITOS2). No external tools.

Run: python3 -m pytest -q tests/test_mito_rescue.py
"""
import random
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_rescue as mr  # noqa: E402

rnd = random.Random(11)
SEQ = Seq("".join(rnd.choice("ACGT") for _ in range(3000)))


def mitoz_rec(feats):
    r = SeqRecord(SEQ, id="S1", description="")
    r.annotations.update(molecule_type="DNA", topology="circular")
    for kind, s, e, strand, gene in feats:
        loc = FeatureLocation(s, e, strand)
        r.features += [SeqFeature(loc, type="gene", qualifiers={"gene": [gene]}),
                       SeqFeature(loc, type=kind, qualifiers={"gene": [gene]})]
    return r


def bed(tmp_path, lines):
    (tmp_path / "result.bed").write_text("".join(f"x\t{s}\t{e}\t{n}\t1\t{st}\n" for s, e, n, st in lines))
    return mr.mitos_to_record(tmp_path, SEQ, "x")


def test_missing_cds_added_mitoz_style(tmp_path):
    """2026-10-20 dev benchmark: atp8 missed by MitoZ (NC_038192.1) but found by MITOS2."""
    r = mitoz_rec([("CDS", 0, 1500, 1, "COX1")])
    m = bed(tmp_path, [(0, 1500, "cox1", "+"), (1600, 1756, "atp8", "-")])
    assert mr.rescue(r, m, 5) == ["ATP8"]
    cds = [f for f in r.features if f.type == "CDS" and f.qualifiers["gene"] == ["ATP8"]]
    assert len(cds) == 1 and cds[0].location.strand == -1
    q = cds[0].qualifiers
    assert q["product"] == ["ATP synthase F0 subunit 8"] and q["transl_table"] == ["5"] and q["codon_start"] == ["1"]
    assert q["translation"][0] == str(cds[0].extract(SEQ).translate(table=5)).rstrip("*")
    assert "MITOS2" in q["note"][0]
    assert any(f.type == "gene" and f.qualifiers["gene"] == ["ATP8"] for f in r.features)


def test_overlapping_prediction_not_added(tmp_path):
    """A MITOS2 gene lying mostly on an existing MitoZ feature is a naming disagreement, not a missed gene."""
    r = mitoz_rec([("CDS", 0, 1500, 1, "COX1")])
    assert mr.rescue(r, bed(tmp_path, [(100, 1400, "nad6", "+")]), 5) == []


def test_rrna_longest_piece_over_minimum(tmp_path):
    """MITOS2 often splits a divergent rRNA into short pieces: only the longest piece >= 300 nt is added."""
    r = mitoz_rec([("CDS", 0, 1000, 1, "COX1")])
    m = bed(tmp_path, [(1100, 1220, "rrnL", "+"), (1300, 1950, "rrnL_1", "+"), (2100, 2350, "rrnS", "+")])
    assert mr.rescue(r, m, 5) == ["l-rRNA"]
    rr = [f for f in r.features if f.type == "rRNA"]
    assert [(int(f.location.start), int(f.location.end)) for f in rr] == [(1300, 1950)]
    assert rr[0].qualifiers["product"] == ["16S ribosomal RNA"]


def test_present_gene_not_duplicated(tmp_path):
    r = mitoz_rec([("rRNA", 0, 800, 1, "s-rRNA")])
    assert mr.rescue(r, bed(tmp_path, [(1000, 1800, "rrnS", "+")]), 5) == []


def test_rescue_dir_in_place(tmp_path):
    """Final GenBank files are updated in place; a report row per added gene."""
    d = tmp_path / "final"
    d.mkdir()
    SeqIO.write(mitoz_rec([("CDS", 0, 1500, 1, "COX1")]), str(d / "S1.gbf"), "genbank")

    def fake_mitos(fa, code, out, linear):
        (out / "result.bed").write_text("x\t1600\t1756\tatp8\t1\t+\n")
        return out

    rows = mr.rescue_dir(d, 5, run=fake_mitos)
    assert [(x["record"], x["gene"]) for x in rows] == [("S1", "ATP8")]
    genes = [f.qualifiers["gene"][0] for f in SeqIO.read(str(d / "S1.gbf"), "genbank").features if f.type == "CDS"]
    assert genes == ["COX1", "ATP8"]


def test_complete_record_skips_mitos2(tmp_path):
    """MITOS2 protein annotation is slow (208 s on NC_038192.1): records with all 13 PCGs + 2 rRNAs are not sent."""
    d = tmp_path / "final"
    d.mkdir()
    names = [v[0] for v in mr.MITOZ_NAMES.values()]
    feats = [("rRNA" if "rRNA" in n else "CDS", 100 * i, 100 * i + 90, 1, n) for i, n in enumerate(names)]
    SeqIO.write(mitoz_rec(feats), str(d / "S1.gbf"), "genbank")
    assert mr.rescue_dir(d, 5, run=lambda *a: (_ for _ in ()).throw(AssertionError("MITOS2 must not run"))) == []
