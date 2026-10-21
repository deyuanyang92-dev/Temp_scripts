"""Tests for mito_rescue.py (genes MitoZ missed, added from MITOS2). No external tools.

Run: python3 -m pytest -q tests/test_mito_rescue.py
"""
import random
import sys
import warnings
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_rescue as mr  # noqa: E402

rnd = random.Random(11)
SEQ = Seq("".join(rnd.choice("ACGT") for _ in range(5000)))


def mitoz_rec(feats, seq=SEQ):
    r = SeqRecord(seq, id="S1", description="")
    r.annotations.update(molecule_type="DNA", topology="circular")
    for kind, s, e, strand, gene in feats:
        loc = FeatureLocation(s, e, strand)
        r.features += [SeqFeature(loc, type="gene", qualifiers={"gene": [gene]}),
                       SeqFeature(loc, type=kind, qualifiers={"gene": [gene]})]
    return r


def bed(tmp_path, lines, seq=SEQ):
    (tmp_path / "result.bed").write_text("".join(f"x\t{s}\t{e}\t{n}\t1\t{st}\n" for s, e, n, st in lines))
    return mr.mitos_to_record(tmp_path, seq, "x")


def test_missing_cds_added_mitoz_style(tmp_path):
    """2026-10-20 dev benchmark: atp8 missed by MitoZ (NC_038192.1) but found by MITOS2.
    (The atp8 region is a clean ORF: a random sequence would carry internal stops and get no /translation.)"""
    orf = Seq("ATG" + "GCT" * 50 + "TAA")                     # 156 nt, minus strand
    seq = SEQ[:1600] + orf.reverse_complement() + SEQ[1756:]
    r = mitoz_rec([("CDS", 0, 1500, 1, "COX1")], seq=seq)
    m = bed(tmp_path, [(0, 1500, "cox1", "+"), (1600, 1756, "atp8", "-")], seq=seq)
    assert mr.rescue(r, m, 5) == ["ATP8"]
    cds = [f for f in r.features if f.type == "CDS" and f.qualifiers["gene"] == ["ATP8"]]
    assert len(cds) == 1 and cds[0].location.strand == -1
    q = cds[0].qualifiers
    assert q["product"] == ["ATP synthase F0 subunit 8"] and q["transl_table"] == ["5"] and q["codon_start"] == ["1"]
    assert q["translation"][0] == "M" + "A" * 50 == str(cds[0].extract(seq).translate(table=5)).rstrip("*")
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


def test_genes_inside_intron_spanning_cds_added(tmp_path):
    """2026-10-21 dev benchmark NC_006902.1 (Metridium-like anthozoan): MitoZ called nad5 as one CDS over
    its ~11 kb group I intron, which carries 10 other genes; all were skipped as 'lying on an existing CDS'.
    A CDS that is much longer than its gene can be (> 1.5 x MAX_LEN) and fully contains >= 2 MITOS2 genes
    of other names does not block them."""
    r = mitoz_rec([("CDS", 100, 4900, 1, "ND5")])
    m = bed(tmp_path, [(100, 400, "nad5", "+"), (500, 1400, "nad1", "+"), (1500, 2600, "cob", "+"),
                       (4600, 4900, "nad5", "+")])
    assert sorted(mr.rescue(r, m, 4)) == ["CYTB", "ND1"]


def test_short_container_never_released(tmp_path):
    """Review C finding 2: ND5 0..1800 (shorter than 1.5 x MAX_LEN nad5) is a correct call; MITOS2 nad4 / nad2
    fragments inside it (homologues) must not be added as two bogus CDS inside a good gene."""
    r = mitoz_rec([("CDS", 0, 1800, 1, "ND5")])
    m = bed(tmp_path, [(100, 700, "nad4", "+"), (800, 1500, "nad2", "+")])
    assert mr.rescue(r, m, 5) == []
    r = mitoz_rec([("CDS", 0, 1800, 1, "ND5")])
    m = bed(tmp_path, [(100, 700, "nad6", "-"), (800, 1100, "nad3", "-")])        # opposite strand
    assert mr.rescue(r, m, 5) == []


def test_container_release_counts_only_absent_genes(tmp_path):
    """Review C finding 2: genes already annotated elsewhere do not count towards releasing a long container."""
    r = mitoz_rec([("CDS", 100, 4900, 1, "ND5"), ("CDS", 4920, 4990, 1, "ND1")])
    m = bed(tmp_path, [(500, 1400, "nad1", "+"), (1500, 2600, "cob", "+")])        # only cob is absent
    assert mr.rescue(r, m, 4) == []


def test_max_len_table_covers_all_genes():
    """Every rescuable gene has an expected maximum length (nt), keyed by MITOS2 name."""
    assert set(mr.MAX_LEN) == set(mr.MITOZ_NAMES)
    assert mr.MAX_LEN["nad5"] == 2000 and mr.MAX_LEN["rrnL"] == 1900 and mr.MAX_LEN["atp8"] == 250


def test_partial_overlap_still_blocks(tmp_path):
    """Only whole contained genes release a CDS: a single overlapping call stays a naming disagreement."""
    r = mitoz_rec([("CDS", 100, 1600, 1, "COX1")])
    assert mr.rescue(r, bed(tmp_path, [(150, 1500, "nad4", "+")]), 5) == []


def test_overlapping_mitos_genes_longest_wins(tmp_path):
    """Review C finding 3: two MITOS2 genes of different names that overlap each other (nad4/nad5 paralogue
    calls) are not both added: longest first, the other lies mostly (> 50 %) on an accepted gene."""
    r = mitoz_rec([("CDS", 0, 500, 1, "COX1")])
    m = bed(tmp_path, [(1000, 1900, "nad4", "+"), (1100, 2200, "nad5", "+")])
    assert mr.rescue(r, m, 5) == ["ND5"]


def test_trna_blocks_overlapping_prediction(tmp_path):
    """Review C finding 3: a MITOS2 gene lying on an annotated tRNA is not added."""
    r = mitoz_rec([("CDS", 0, 500, 1, "COX1"), ("tRNA", 1000, 1070, 1, "trnX")])
    assert mr.rescue(r, bed(tmp_path, [(1000, 1070, "atp8", "+")]), 5) == []


def test_rrna_fragments_not_joined(tmp_path):
    """Review C finding 5: two short rRNA pieces far apart are not one spliced rRNA (sum of pieces >= 300)."""
    r = mitoz_rec([("CDS", 0, 500, 1, "COX1")])
    m = bed(tmp_path, [(100, 260, "rrnL-a", "+"), (3000, 3160, "rrnL-b", "+")])
    assert [f.type for f in m.features].count("rRNA") == 2          # kept as separate features
    assert mr.rescue(r, m, 5) == []


def test_rrna_fragment_longest_piece_added(tmp_path):
    """rRNA fragments stay separate; the longest piece >= MIN_RRNA is the one added (as the docstring says)."""
    r = mitoz_rec([("CDS", 0, 500, 1, "COX1")])
    m = bed(tmp_path, [(1000, 1200, "rrnL-a", "+"), (1300, 1700, "rrnL-b", "+")])
    assert mr.rescue(r, m, 5) == ["l-rRNA"]
    rr = [f for f in r.features if f.type == "rRNA"]
    assert [(int(f.location.start), int(f.location.end)) for f in rr] == [(1300, 1700)]


def test_cds_fragments_joined_only_when_close_and_same_strand(tmp_path):
    """Review C finding 5: CDS fragments (-a/-b) are joined when on the same strand with a gap <= 3000 nt."""
    near = bed(tmp_path, [(1000, 1300, "cox2-a", "+"), (2000, 2300, "cox2-b", "+")])
    cds = [f for f in near.features if f.type == "CDS"]
    assert len(cds) == 1 and [(int(p.start), int(p.end)) for p in cds[0].location.parts] == [(1000, 1300), (2000, 2300)]
    far = bed(tmp_path, [(100, 400, "cox2-a", "+"), (3600, 3900, "cox2-b", "+")])         # gap 3200
    assert len([f for f in far.features if f.type == "CDS"]) == 2
    opp = bed(tmp_path, [(1000, 1300, "cox2-a", "+"), (2000, 2300, "cox2-b", "-")])
    assert len([f for f in opp.features if f.type == "CDS"]) == 2


def _tx_case(tmp_path, cds_nt, code, tail=""):
    seq = Seq("AAAAAAAAAA" + cds_nt + tail + "C" * 20)
    r = mitoz_rec([("CDS", 0, 5, 1, "COX1")], seq=seq)
    m = bed(tmp_path, [(10, 10 + len(cds_nt) + len(tail), "nad6", "+")], seq=seq)
    with warnings.catch_warnings():
        warnings.simplefilter("error")                  # BiopythonWarning (partial codon) must not occur
        assert mr.rescue(r, m, code) == ["ND6"]
    return [f for f in r.features if f.type == "CDS" and f.qualifiers["gene"] == ["ND6"]][0]


def test_translation_trimmed_alt_start_is_m(tmp_path):
    """Review C finding 4: length not a multiple of 3 -> trimmed, no warning; table 5 start ATT -> M;
    trailing stop removed."""
    f = _tx_case(tmp_path, "ATT" + "GCC" * 20 + "TAA", 5, tail="CC")
    assert f.qualifiers["translation"] == ["M" + "A" * 20]


def test_translation_internal_stop_omitted(tmp_path):
    """Review C finding 4: an internal stop means the CDS is doubtful: no /translation is written."""
    f = _tx_case(tmp_path, "ATG" + "GCC" * 5 + "TAA" + "GCC" * 5 + "TAA", 5)
    assert "translation" not in f.qualifiers and f.qualifiers["gene"] == ["ND6"]


def test_translation_non_start_first_codon_kept(tmp_path):
    """A first codon that is not a start codon of the table is translated as is."""
    f = _tx_case(tmp_path, "GCC" + "GCC" * 10 + "TAA", 5)
    assert f.qualifiers["translation"] == ["A" * 11]


def test_note_is_neutral_and_hint_goes_to_report(tmp_path):
    """Review C finding 15: no tool name / to-do in the feature note; the hint is a report column."""
    d = tmp_path / "final"
    d.mkdir()
    SeqIO.write(mitoz_rec([("CDS", 0, 1500, 1, "COX1")]), str(d / "S1.gbf"), "genbank")

    def fake_mitos(fa, code, out, linear):
        (out / "result.bed").write_text("x\t1600\t1756\tatp8\t1\t+\n")
        return out

    rows = mr.rescue_dir(d, 5, run=fake_mitos)
    assert "boundar" in rows[0]["hint"]
    notes = [n for f in SeqIO.read(str(d / "S1.gbf"), "genbank").features if f.type == "CDS"
             for n in f.qualifiers.get("note", [])]
    assert notes == ["annotated by similarity to a MITOS2 prediction (not found by the primary annotation)"]
    assert mr.NOTE == notes[0]


def test_report_has_hint_column(tmp_path, monkeypatch):
    d = tmp_path / "final"
    d.mkdir()
    SeqIO.write(mitoz_rec([("CDS", 0, 1500, 1, "COX1")]), str(d / "S1.gbf"), "genbank")

    def fake_run_mitos(fa, code, out, linear, **kw):
        (out / "result.bed").write_text("x\t1600\t1756\tatp8\t1\t+\n")
        return out

    monkeypatch.setattr(mr, "run_mitos", fake_run_mitos)
    assert mr.main(["-i", str(d), "--genetic_code", "5"]) == 0
    lines = (d / "rescue_report.tsv").read_text().splitlines()
    head = lines[0].split("\t")
    assert head[:5] == ["file", "record", "gene", "location", "status"] and "hint" in head
    row = dict(zip(head, lines[1].split("\t")))
    assert row["gene"] == "ATP8" and row["status"] == "added" and "boundar" in row["hint"]
