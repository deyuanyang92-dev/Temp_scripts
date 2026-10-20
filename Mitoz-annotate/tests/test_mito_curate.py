"""Tests for mito_curate.py on real MitoZ output (Aphrodita KY753833; Lepidonotus clava OW387151).

Rules: Fourdrilis et al. 2018 criteria as summarised by Ghiselli et al. 2021 (doi:10.1098/rstb.2020.0159).
Run: python3 -m pytest -q tests/test_mito_curate.py
"""
import csv
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
SCRIPT = HERE.parent / "mito_curate.py"


def curate(tmp_path, fixture, *extra):
    (tmp_path / "in").mkdir(exist_ok=True)
    (tmp_path / "in" / fixture).write_bytes((DATA / fixture).read_bytes())
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(tmp_path / "in"), "-o",
                    str(tmp_path / "out"), "--genetic_code", "5", *extra], check=True, capture_output=True)
    rec = SeqIO.read(str(tmp_path / "out" / fixture), "genbank")
    rows = list(csv.DictReader((tmp_path / "out" / "curation_report.tsv").open(), delimiter="\t"))
    return rec, rows


def feat(rec, ftype, start1):
    return next(f for f in rec.features if f.type == ftype and int(f.location.start) + 1 == start1)


def test_rrna_trimmed_not_overlapping_trnas(tmp_path):
    rec, _ = curate(tmp_path, "KY753833_mitoz.gbf")
    rrna = sorted((int(f.location.start) + 1, int(f.location.end), f.qualifiers["gene"][0], f.qualifiers["product"][0])
                  for f in rec.features if f.type == "rRNA")
    assert rrna == [(9833, 10643, "rrnS", "12S ribosomal RNA"), (10708, 11997, "rrnL", "16S ribosomal RNA")]
    trnas = [(int(f.location.start) + 1, int(f.location.end)) for f in rec.features if f.type == "tRNA"]
    for s, e, *_ in rrna:
        assert not any(ts <= e and te >= s for ts, te in trnas)          # criterion vi: no overlap


def test_rrna_fill_extends_into_free_space_only(tmp_path):
    rec, _ = curate(tmp_path, "KY753833_mitoz.gbf", "--rrna_bounds", "fill")
    rrna = {f.qualifiers["gene"][0]: (int(f.location.start) + 1, int(f.location.end)) for f in rec.features if f.type == "rRNA"}
    assert rrna["rrnS"] == (9833, 10643) and rrna["rrnL"] == (10708, 11997)  # already flank-to-flank here


def test_cox1_abbreviated_stop_before_trna(tmp_path):
    """MitoZ cox1 1..1543 has no stop; the first complete stop lies inside trnN (1557..): the gene
    ends with TA right before trnN (NCBI KY753833 likewise ends COX1 with TA before tRNA-Asn)."""
    rec, rows = curate(tmp_path, "KY753833_mitoz.gbf")
    cox1 = feat(rec, "CDS", 1)
    assert int(cox1.location.end) == 1556
    assert cox1.qualifiers["transl_except"] == ["(pos:1555..1556,aa:TERM)"]
    assert any(r["rule"] == "P3" and r["gene"].upper() == "COX1" for r in rows)


def test_different_frame_overlap_keeps_start(tmp_path):
    """7994..8281 and 8275..9642 overlap in different frames (criterion ii allows it): start kept."""
    rec, _ = curate(tmp_path, "KY753833_mitoz.gbf")
    assert int(feat(rec, "CDS", 8275).location.start) + 1 == 8275


def test_internal_stops_reported_not_changed(tmp_path):
    before = SeqIO.read(str(DATA / "OW387151_mitoz.gbf"), "genbank")
    rec, rows = curate(tmp_path, "OW387151_mitoz.gbf")
    bad = [r for r in rows if r["rule"] == "P4" and "internal stop" in r["reason"]]
    assert bad and bad[0]["gene"].upper() == "COX1"
    loc = {str(f.location) for f in before.features if f.type == "CDS"}
    assert any(str(f.location) in loc for f in rec.features if f.type == "CDS" and f.qualifiers["gene"][0].upper() == "COX1")


def test_ncbi_validation_before_after(tmp_path):
    """Built-in NCBI table2asn check before vs after curation: errors fixed, none introduced."""
    import shutil
    t2a = shutil.which("table2asn") or "/root/miniconda3/envs/mitointron/bin/table2asn"
    if not Path(t2a).exists():
        import pytest
        pytest.skip("table2asn not available")
    curate(tmp_path, "KY753833_mitoz.gbf", "--table2asn", t2a)
    rows = list(csv.DictReader((tmp_path / "out" / "validation_compare.tsv").open(), delimiter="\t"))
    err = [r for r in rows if r["severity"] == "ERROR"]
    assert any(r["status"] == "fixed" for r in err)
    assert not any(r["status"] == "new" for r in err)


def _spacer_record(tmp_path, spacer):
    """+ strand gene ATG(GCC)x30 + T, then `spacer`, then a 66-nt tRNA whose frame holds the first TAA."""
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    gene = "ATG" + "GCC" * 30 + "T"                                  # 94 nt, abbreviated stop T
    k = -(len(gene) + len(spacer) + 30) % 3                          # put the TAA in the gene's frame
    trna = ("GGCCGC" * 5 + "G" * k + "TAA" + "GCGGCC" * 5 + "CCC")[:66]
    pad = "CCCGGG" * 20
    seq = gene + spacer + trna + pad
    t0 = len(gene) + len(spacer)
    frame_stop = next(i for i in range(0, len(seq) - 2, 3) if seq[i:i + 3] in ("TAA", "TAG"))
    rec = SeqRecord(Seq(seq), id="X", description="")
    rec.annotations.update(molecule_type="DNA", topology="linear")
    cds = FeatureLocation(0, frame_stop + 3, 1)                      # MitoZ: greedy stop inside the tRNA
    rec.features = [SeqFeature(FeatureLocation(0, len(seq), 1), type="source"),
                    SeqFeature(cds, type="gene", qualifiers={"gene": ["COX3"]}),
                    SeqFeature(cds, type="CDS", qualifiers={"gene": ["COX3"], "codon_start": ["1"], "transl_table": ["5"]}),
                    SeqFeature(FeatureLocation(t0, t0 + len(trna), 1), type="tRNA",
                               qualifiers={"gene": ["trnS(uga)"], "product": ["tRNA-Ser"]})]
    d = tmp_path / "in"
    d.mkdir()
    SeqIO.write(rec, str(d / "x.gbf"), "genbank")
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(d), "-o", str(tmp_path / "out"),
                    "--genetic_code", "5", "--no_validate"], check=True, capture_output=True)
    out = SeqIO.read(str(tmp_path / "out" / "x.gbf"), "genbank")
    return next(f for f in out.features if f.type == "CDS")


def test_abbreviated_stop_with_short_spacer(tmp_path):
    """2026-10-20 dev benchmark: NCBI ends incomplete-stop genes 0 (143x), -1 (9x), 1-3 nt (10x) before the next
    same-strand tRNA; P3 only accepted 0, so e.g. NC_018362.1 nad2 (T, 1 nt, trnE) kept MitoZ's stop inside the tRNA."""
    cds = _spacer_record(tmp_path, "CC")
    assert int(cds.location.end) == 94 and cds.qualifiers["transl_except"] == ["(pos:94,aa:TERM)"]


def test_abbreviated_stop_no_spacer_unchanged(tmp_path):
    assert int(_spacer_record(tmp_path, "").location.end) == 94


def test_rrna_never_contains_a_trna(tmp_path):
    """2026-10-21 table2asn on dev outputs: 16x SEQ_FEAT.BadRRNAcomponentOverlapTRNA, e.g. NC_028229.1 rrnS 3368..4148
    containing trnE 3373..3429 - R1 only trimmed neighbours overlapping an rRNA end. The rRNA keeps the longer side."""
    import random
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    rnd = random.Random(2)
    rec = SeqRecord(Seq("".join(rnd.choice("ACGT") for _ in range(2000))), id="X", description="")
    rec.annotations.update(molecule_type="DNA", topology="linear")
    def add(kind, s, e, gene, product):
        loc = FeatureLocation(s - 1, e, 1)
        rec.features += [SeqFeature(loc, type="gene", qualifiers={"gene": [gene]}),
                         SeqFeature(loc, type=kind, qualifiers={"gene": [gene], "product": [product]})]
    rec.features = [SeqFeature(FeatureLocation(0, 2000, 1), type="source")]
    add("tRNA", 100, 160, "trnA(ugc)", "tRNA-Ala")
    add("rRNA", 150, 1000, "s-rRNA", "12S ribosomal RNA")          # overlaps trnA end, contains trnE
    add("tRNA", 170, 230, "trnE(uuc)", "tRNA-Glu")
    add("tRNA", 1001, 1060, "trnS(uga)", "tRNA-Ser")
    d = tmp_path / "in"
    d.mkdir()
    SeqIO.write(rec, str(d / "x.gbf"), "genbank")
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(d), "-o", str(tmp_path / "out"),
                    "--genetic_code", "5", "--no_validate"], check=True, capture_output=True)
    out = SeqIO.read(str(tmp_path / "out" / "x.gbf"), "genbank")
    r = next(f for f in out.features if f.type == "rRNA")
    assert (int(r.location.start) + 1, int(r.location.end)) == (231, 1000)
    assert r.qualifiers["gene"] == ["rrnS"]


def test_ncbi_gene_names(tmp_path):
    """User 2026-10-21: names consistent and NCBI-conform. /gene is the plain symbol (NCBI genome annotation guide:
    tRNA labels go to /note); Leu/Ser by anticodon (L1 CUN, L2 UUR, S1 AGN, S2 UCN); rRNA rrnS/rrnL also when the
    feature is a join (3 dev records kept MitoZ's l-rRNA / s-rRNA); CDS synonyms -> COX1..CYTB, ND1..ND6, ND4L."""
    import random
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    rnd = random.Random(3)
    rec = SeqRecord(Seq("".join(rnd.choice("ACGT") for _ in range(3000))), id="X", description="")
    rec.annotations.update(molecule_type="DNA", topology="circular")
    rec.features = [SeqFeature(FeatureLocation(0, 3000, 1), type="source")]
    def add(kind, loc, gene, **q):
        rec.features += [SeqFeature(loc, type="gene", qualifiers={"gene": [gene]}),
                         SeqFeature(loc, type=kind, qualifiers={"gene": [gene], **{k: [v] for k, v in q.items()}})]
    add("tRNA", FeatureLocation(100, 165, 1), "trnL(uaa)", product="tRNA-Leu")
    add("tRNA", FeatureLocation(200, 265, -1), "trnS(gcu)", product="tRNA-Ser")
    add("tRNA", FeatureLocation(300, 365, 1), "trnM(cau)", product="tRNA-Met")
    add("rRNA", CompoundLocation([FeatureLocation(2500, 3000, 1), FeatureLocation(0, 90, 1)]), "l-rRNA",
        product="16S ribosomal RNA")
    add("CDS", FeatureLocation(400, 1000, 1), "cob", product="Cytb", codon_start="1", transl_table="5")
    d = tmp_path / "in"
    d.mkdir()
    SeqIO.write(rec, str(d / "x.gbf"), "genbank")
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(d), "-o", str(tmp_path / "out"),
                    "--genetic_code", "5", "--no_validate"], check=True, capture_output=True)
    out = SeqIO.read(str(tmp_path / "out" / "x.gbf"), "genbank")
    got = [(f.type, f.qualifiers["gene"][0], f.qualifiers.get("product", [""])[0], f.qualifiers.get("note", [""])[0])
           for f in out.features if f.type in ("tRNA", "rRNA", "CDS")]
    assert ("tRNA", "trnL2", "tRNA-Leu", "anticodon: UAA") in got
    assert ("tRNA", "trnS1", "tRNA-Ser", "anticodon: GCU") in got
    assert ("tRNA", "trnM", "tRNA-Met", "anticodon: CAU") in got
    assert ("rRNA", "rrnL", "16S ribosomal RNA", "") in got
    assert any(x[:3] == ("CDS", "CYTB", "cytochrome b") for x in got)
    genes = sorted(f.qualifiers["gene"][0] for f in out.features if f.type == "gene")
    assert genes == ["CYTB", "rrnL", "trnL2", "trnM", "trnS1"]


def _one_cds(tmp_path, nt, code, gene_loc=None):
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    seq = "CCCGGG" * 10 + nt + "CCCGGG" * 10
    rec = SeqRecord(Seq(seq), id="X", description="")
    rec.annotations.update(molecule_type="DNA", topology="linear")
    loc = FeatureLocation(60, 60 + len(nt), 1)
    rec.features = [SeqFeature(FeatureLocation(0, len(seq), 1), type="source"),
                    SeqFeature(gene_loc or loc, type="gene", qualifiers={"gene": ["ND4L"]}),
                    SeqFeature(loc, type="CDS", qualifiers={"gene": ["ND4L"], "codon_start": ["1"], "transl_table": [str(code)]})]
    d = tmp_path / "in"
    d.mkdir()
    SeqIO.write(rec, str(d / "x.gbf"), "genbank")
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(d), "-o", str(tmp_path / "out"),
                    "--genetic_code", str(code), "--no_validate"], check=True, capture_output=True)
    return SeqIO.read(str(tmp_path / "out" / "x.gbf"), "genbank")


def test_illegal_start_marked_partial(tmp_path):
    """2026-10-21 table2asn: SEQ_FEAT.StartCodon ERROR on echinoderm nad4l (ATT/ATC, code 9 starts are ATG/GTG only);
    NCBI writes the same CDS as <start with note 'start codon not determined' (NC_001627.1, NC_061414.1)."""
    out = _one_cds(tmp_path, "ATT" + "GCC" * 40 + "TAA", 9)
    cds = next(f for f in out.features if f.type == "CDS")
    assert str(cds.location).startswith("[<60") and "start codon not determined" in cds.qualifiers["note"]
    gene = next(f for f in out.features if f.type == "gene")
    assert str(gene.location) == str(cds.location)


def test_gene_feature_matches_cds(tmp_path):
    """table2asn SEQ_FEAT.CDSgeneRange: gene features must span exactly their CDS / RNA."""
    from Bio.SeqFeature import FeatureLocation
    out = _one_cds(tmp_path, "ATG" + "GCC" * 40 + "TAA", 5, gene_loc=FeatureLocation(55, 150, 1))
    cds, gene = (next(f for f in out.features if f.type == t) for t in ("CDS", "gene"))
    assert (int(gene.location.start), int(gene.location.end)) == (int(cds.location.start), int(cds.location.end))


def test_circular_defline_complete():
    """table2asn SEQ_INST.CompleteCircleProblem on every circular record: the .fsa needs [completeness=complete]."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    sys.path.insert(0, str(HERE.parent))
    import ncbi_validate as nv
    r = SeqRecord(Seq("ACGT"), id="X")
    r.annotations.update(organism="Aus bus", topology="circular")
    assert "[completeness=complete]" in nv.fasta_defline(r, 5)
    r.annotations["topology"] = "linear"
    assert "completeness" not in nv.fasta_defline(r, 5)


def test_missing_stop_marked_partial(tmp_path):
    """2026-10-21 dev benchmark: SEQ_FEAT.NoStop ERROR (e.g. NC_044095.1 ND6, NC_085563.1 ND6) when a CDS has neither
    a complete nor an abbreviated (T/TA) stop; written 3'-partial with 'stop codon not determined' like P5."""
    out = _one_cds(tmp_path, "ATG" + "GCC" * 40, 5)
    cds = next(f for f in out.features if f.type == "CDS")
    assert str(cds.location).endswith(">183](+)") and "stop codon not determined" in cds.qualifiers["note"]


def _origin_record(tmp_path, strand):
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    body = "CCCGGG" * 20
    gene = "ATG" + "GCC" * 21 + "TAA"                      # 69 nt; the last 4 nt lie past the origin
    plus = gene[-4:] + body + gene[:-4]
    seq = plus if strand == 1 else str(Seq(plus).reverse_complement())
    L = len(seq)
    loc = FeatureLocation(L - 65, L - 2, 1) if strand == 1 else FeatureLocation(2, 65, -1)
    rec = SeqRecord(Seq(seq), id="X", description="")
    rec.annotations.update(molecule_type="DNA", topology="circular")
    rec.features = [SeqFeature(FeatureLocation(0, L, 1), type="source"),
                    SeqFeature(loc, type="gene", qualifiers={"gene": ["ATP8"]}),
                    SeqFeature(loc, type="CDS", qualifiers={"gene": ["ATP8"], "codon_start": ["1"], "transl_table": ["5"]})]
    d = tmp_path / "in"
    d.mkdir(parents=True)
    SeqIO.write(rec, str(d / "x.gbf"), "genbank")
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(d), "-o", str(tmp_path / "out"),
                    "--genetic_code", "5", "--no_validate"], check=True, capture_output=True)
    return SeqIO.read(str(tmp_path / "out" / "x.gbf"), "genbank"), gene


def test_cds_extended_across_origin(tmp_path):
    """2026-10-21 dev benchmark NC_006902.1 / NC_022829.1 (ATP8 overlaps COX1, cut by the cox1-first rotation),
    NC_025335.1 (ND4L): CDS stops 2 nt before the end of a circular record; its stop lies past the origin."""
    for strand in (1, -1):
        out, gene = _origin_record(tmp_path / str(strand), strand)
        cds = next(f for f in out.features if f.type == "CDS")
        assert len(cds.location.parts) == 2 and str(cds.extract(out.seq)) == gene, strand
        g = next(f for f in out.features if f.type == "gene")
        assert str(g.location) == str(cds.location)


def test_tbl_minus_strand_join_in_biological_order():
    """2026-10-21: NC_035061.1 trnY complement(join(16132..16186,1..8)) was written to the .tbl as 16186-16132
    then 8-1; a feature table lists intervals in transcript order (8-1, then 16186-16132). Biopython already
    keeps minus-strand parts in transcript order."""
    import io
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    sys.path.insert(0, str(HERE.parent))
    import ncbi_validate as nv
    r = SeqRecord(Seq("A" * 100), id="X")
    r.annotations["molecule_type"] = "DNA"
    r.features = [SeqFeature(CompoundLocation([FeatureLocation(0, 10, -1), FeatureLocation(90, 100, -1)]), type="tRNA")]
    h = io.StringIO()
    SeqIO.write(r, h, "genbank")
    assert "complement(join(91..100,1..10))" in h.getvalue()
    assert nv.feature_lines(r)[1:3] == ["10\t1\ttRNA", "100\t91"]
