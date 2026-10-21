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


# ── in-process helpers (review C fixes) ──
def _mc():
    sys.path.insert(0, str(HERE.parent))
    import mito_curate
    return mito_curate


def _args(**kw):
    import argparse
    d = dict(genetic_code=5, rrna_bounds="fill", max_fill=300, max_start_shift=30, aggressive_start=False)
    d.update(kw)
    return argparse.Namespace(**d)


def _mkrec(seq, feats, topology="linear"):
    """feats: (type, location, gene[, extra qualifiers]); a matching gene feature is added for CDS / RNA."""
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(Seq(seq), id="X", description="")
    rec.annotations.update(molecule_type="DNA", topology=topology)
    rec.features = [SeqFeature(FeatureLocation(0, len(seq), 1), type="source")]
    for ftype, loc, gene, *extra in feats:
        q = {"gene": [gene], **(extra[0] if extra else {})}
        if ftype == "CDS":
            q.setdefault("codon_start", ["1"])
            q.setdefault("transl_table", ["5"])
        rec.features += [SeqFeature(loc, type="gene", qualifiers={"gene": [gene]}), SeqFeature(loc, type=ftype, qualifiers=q)]
    return rec


def _curate_rec(rec, **kw):
    a = _args(**kw)
    rows = _mc().Curator(rec, a.genetic_code, a).run()
    return rows


def _by(rec, ftype, gene):
    return next(f for f in rec.features if f.type == ftype and f.qualifiers.get("gene", [""])[0] == gene)


def test_wrapped_trna_does_not_block_p1():
    """Review C finding 6: span() is the hull of a compound location, so an origin-spanning tRNA
    join(L-29..L,1..35) counted as a preceding feature ending at L and disabled P1 for every same-strand CDS."""
    from Bio.SeqFeature import CompoundLocation, FeatureLocation
    L = 1200
    seq = "GCC" * 100 + "ATG" + "GCC" + "GCC" * 40 + "TAA" + "GCC" * 257        # ATG 301..303, annotated start 307 = GCC
    assert len(seq) == L
    cds = FeatureLocation(306, 429, 1)
    trna = CompoundLocation([FeatureLocation(L - 30, L, 1), FeatureLocation(0, 35, 1)])
    rec = _mkrec(seq, [("tRNA", trna, "trnY", {"product": ["tRNA-Tyr"]}), ("CDS", cds, "ND1")], "circular")
    rows = _curate_rec(rec)
    assert any(r["rule"] == "P1" for r in rows), rows
    assert int(_by(rec, "CDS", "ND1").location.start) + 1 == 301


def _rrna_rec(feats, n=3000):
    import random
    rnd = random.Random(11)
    seq = "".join(rnd.choice("ACGT") for _ in range(n))
    from Bio.SeqFeature import FeatureLocation
    out = []
    for kind, s, e, gene in feats:
        out.append((kind, FeatureLocation(s - 1, e, 1), gene, {"product": [{"tRNA": "tRNA-Gln", "rRNA": "16S ribosomal RNA"}[kind]]}))
    return _mkrec(seq, out)


def test_rrna_not_trimmed_for_inside_gene_when_loss_large():
    """Review C finding 7: a spurious tRNA call inside a correctly annotated 1.3-kb l-rRNA cut the rRNA to 830 nt
    (36 % lost). Trim only when the discarded side is <= 10 % of the rRNA; else keep it and log R3."""
    rec = _rrna_rec([("rRNA", 1001, 2300, "l-rRNA"), ("tRNA", 1401, 1470, "trnQ(uug)")])
    rows = _curate_rec(rec)
    r = _by(rec, "rRNA", "rrnL")
    assert (int(r.location.start) + 1, int(r.location.end)) == (1001, 2300)
    r3 = [x for x in rows if x["rule"] == "R3" and "inside rRNA" in x["reason"]]
    assert len(r3) == 1 and "trnQ" in r3[0]["reason"] and "not trimmed" in r3[0]["reason"], rows
    assert not any(x["rule"] == "R1" for x in rows)


def test_rrna_trimmed_for_inside_gene_near_the_end():
    rec = _rrna_rec([("rRNA", 1001, 2300, "l-rRNA"), ("tRNA", 2201, 2270, "trnQ(uug)")])
    rows = _curate_rec(rec)
    r = _by(rec, "rRNA", "rrnL")
    assert (int(r.location.start) + 1, int(r.location.end)) == (1001, 2200)
    assert any(x["rule"] == "R1" for x in rows)


def test_translation_kept_for_every_cds(tmp_path):
    """Review C finding 9: pcg() popped /translation from every simple CDS and main() now always rewrites the file, so
    KY753833 went from 13 /translation qualifiers to 0 (changed or not). Recomputed from the final location."""
    before = SeqIO.read(str(DATA / "KY753833_mitoz.gbf"), "genbank")
    rec, _ = curate(tmp_path, "KY753833_mitoz.gbf")
    n_in = sum("translation" in f.qualifiers for f in before.features if f.type == "CDS")
    cds = [f for f in rec.features if f.type == "CDS"]
    assert n_in == 13 and sum("translation" in f.qualifiers for f in cds) == 13
    for f in cds:
        prot = f.qualifiers["translation"][0]
        assert prot.startswith("M") and "*" not in prot, f.qualifiers["gene"]
        assert prot == str(f.extract(rec.seq)[: len(f.extract(rec.seq)) // 3 * 3].translate(table=5)).rstrip("*"), \
            f.qualifiers["gene"]


def test_translation_alternative_start_is_methionine():
    """NCBI convention: a (non-partial) CDS beginning with a start codon of the table (ATT under code 5) is
    translated with M; checked on a + CDS, a - CDS and a join."""
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation
    gene = "ATT" + "GCC" * 40 + "TAA"                                       # 126 nt
    pad = "CCCGGG" * 10
    ex1, ex2 = "ATT" + "GCC" * 20, "GCC" * 20 + "TAA"
    seq = pad + gene + pad + str(Seq(gene).reverse_complement()) + pad + ex1 + "A" * 20 + ex2 + pad
    j = 60 + 126 + 60 + 126 + 60                                             # exon 1 starts here (0-based)
    rec = _mkrec(seq, [("CDS", FeatureLocation(60, 186, 1), "ND1"),
                       ("CDS", FeatureLocation(246, 372, -1), "ND2"),
                       ("CDS", CompoundLocation([FeatureLocation(j, j + 63, 1), FeatureLocation(j + 83, j + 146, 1)]), "ND3")])
    _curate_rec(rec)
    for gene_name in ("ND1", "ND2", "ND3"):
        assert _by(rec, "CDS", gene_name).qualifiers["translation"] == ["M" + "A" * 40], gene_name


def test_translation_not_written_with_internal_stop_and_partial_start_plain():
    """An internal stop remains -> no /translation (P4 reports it); a 5'-partial CDS keeps the plain translation (I)."""
    from Bio.SeqFeature import BeforePosition, FeatureLocation
    nt = "ATG" + "GCC" * 10 + "TAA" + "GCC" * 100 + "TAA"
    rec = _mkrec("CCCGGG" * 10 + nt + "CCCGGG" * 10, [("CDS", FeatureLocation(60, 60 + len(nt), 1), "ND1")])
    rows = _curate_rec(rec)
    assert "translation" not in _by(rec, "CDS", "ND1").qualifiers
    assert any(r["rule"] == "P4" and "internal stop" in r["reason"] for r in rows)
    nt = "ATT" + "GCC" * 20 + "TAA"
    rec = _mkrec("CCCGGG" * 10 + nt + "CCCGGG" * 10, [("CDS", FeatureLocation(BeforePosition(60), 60 + len(nt), 1), "ND1")])
    _curate_rec(rec)
    assert _by(rec, "CDS", "ND1").qualifiers["translation"] == ["I" + "A" * 20]


def test_p6_stop_tested_in_frame_only():
    """Review C finding 10: an out-of-frame TAA at the end of a CDS whose length is not a multiple of 3 suppressed P6
    (table2asn then gives SEQ_FEAT.NoStop). Stop test only when length % 3 == 0; otherwise T / TA tail or no stop."""
    from Bio.SeqFeature import FeatureLocation
    nt = "ATG" + "GCC" * 40 + "ATAA"                                       # 127 nt, no in-frame stop
    rec = _mkrec("CCCGGG" * 10 + nt + "CCCGGG" * 10, [("CDS", FeatureLocation(60, 60 + len(nt), 1), "ND1")])
    rows = _curate_rec(rec)
    assert any(r["rule"] == "P6" for r in rows), rows
    assert any(r["rule"] == "P4" and "not a multiple of 3" in r["reason"] for r in rows), rows
    assert "stop codon not determined" in _by(rec, "CDS", "ND1").qualifiers["note"]


def _cut_cds_record(stop_at, first_feature=True):
    """Circular record; + strand ND4L = ATG + GCC*40 ends at L-1 (no stop); its first in-frame stop lies past the origin,
    TAA at 0-based s[stop_at..] (stop_at % 3 == 2). First feature COX1 covers position 1 (cox1-first rotation)."""
    from Bio.SeqFeature import FeatureLocation
    first = list(("GCC" * 134)[:400])
    first[stop_at:stop_at + 3] = "TAA"
    nd4l = "ATG" + "GCC" * 40
    seq = "".join(first) + "CCCGGG" * 30 + nd4l + "G"
    feats = [("CDS", FeatureLocation(580, 703, 1), "ND4L")]
    if first_feature:
        feats.append(("tRNA", FeatureLocation(0, 400, 1), "trnY", {"product": ["tRNA-Tyr"]}))
    rec = _mkrec(seq, feats, "circular")
    return rec, _curate_rec(rec)


def test_origin_extension_capped_by_feature_at_position_1():
    """Review C finding 11: O1 extended a CDS up to 300 nt past the origin, i.e. up to 300 nt into the feature at
    position 1 (cox1-first rotation); genuine mitochondrial overlaps are <= ~50 nt (dev cases 13-19 nt)."""
    rec, rows = _cut_cds_record(14)                               # stop 17 nt into the first feature: genuine overlap
    assert any(r["rule"] == "O1" for r in rows)
    assert len(_by(rec, "CDS", "ND4L").location.parts) == 2
    rec, rows = _cut_cds_record(149)                              # stop 152 nt into the first feature: not O1
    assert not any(r["rule"] == "O1" for r in rows), rows
    nd4l = _by(rec, "CDS", "ND4L")
    assert len(nd4l.location.parts) == 1 and str(nd4l.location).endswith(">703](+)")
    assert any(r["rule"] == "P6" for r in rows)
    rec, rows = _cut_cds_record(149, first_feature=False)         # nothing at position 1: the 300 nt limit applies
    assert any(r["rule"] == "O1" for r in rows)


def test_cds_extension_does_not_overlap_following_rrna():
    """Review C finding 12: R1 ran before the PCG loop, so P3 later extended a CDS (no stop at 1453, first stop at 1531)
    into the rRNA that R1 had just cleared (78-nt CDS/rRNA overlap). The PCG loop now runs first."""
    from Bio.SeqFeature import FeatureLocation
    seq = ("CCCGGG" * 167)[:1000] + "ATG" + "GCC" * 150 + "GCC" * 25 + "TAA" + ("CCCGGG" * 300)[:1469]
    assert len(seq) == 3000 and seq[1528:1531] == "TAA"
    rec = _mkrec(seq, [("CDS", FeatureLocation(1000, 1453, 1), "COX1"),
                       ("rRNA", FeatureLocation(1453, 2500, 1), "s-rRNA", {"product": ["12S ribosomal RNA"]})])
    rows = _curate_rec(rec)
    cds, rrna = _by(rec, "CDS", "COX1"), _by(rec, "rRNA", "rrnS")
    assert int(cds.location.end) == 1531 and any(r["rule"] == "P3" for r in rows)
    assert int(rrna.location.start) >= int(cds.location.end)


def test_failed_curation_leaves_record_unchanged(tmp_path, monkeypatch):
    """Review C finding 14: Curator.run() mutates the record in place, so after an exception the row said 'left
    unchanged' but renames / trims done before the failure were still written."""
    mc = _mc()
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "k.gbf").write_bytes((DATA / "KY753833_mitoz.gbf").read_bytes())

    def boom(self):
        self.names(next(f for f in self.rec.features if f.type == "CDS"))
        self.rec.features[1].qualifiers["gene"] = ["MUTATED"]
        self.rec.features.pop()
        raise RuntimeError("boom")

    monkeypatch.setattr(mc.Curator, "run", boom)
    assert mc.main(["-i", str(tmp_path / "in"), "-o", str(tmp_path / "out"), "--genetic_code", "5", "--no_validate"]) == 0
    before, after = (SeqIO.read(str(tmp_path / d / "k.gbf"), "genbank") for d in ("in", "out"))

    def sig(r):
        return [(f.type, str(f.location), f.qualifiers.get("gene")) for f in r.features]
    assert sig(after) == sig(before)
    rows = list(csv.DictReader((tmp_path / "out" / "curation_report.tsv").open(), delimiter="\t"))
    assert any(r["rule"] == "error" and "left unchanged" in r["reason"] for r in rows)


def test_origin_spanning_feature_sorts_by_its_larger_coordinate_part():
    """Review C finding 20: sort key int(location.start) put join(L-29..L,1..35) (start 0) before every feature."""
    from Bio.SeqFeature import CompoundLocation, FeatureLocation
    seq = "CCCGGG" * 100 + "ATG" + "GCC" * 20 + "TAA" + "CCCGGG" * 100
    L = len(seq)
    trna = CompoundLocation([FeatureLocation(L - 30, L, 1), FeatureLocation(0, 35, 1)])
    nt = 3 + 60 + 3
    rec = _mkrec(seq, [("tRNA", trna, "trnY", {"product": ["tRNA-Tyr"]}), ("CDS", FeatureLocation(600, 600 + nt, 1), "ND1")],
                 "circular")
    _curate_rec(rec)
    order = [(f.type, f.qualifiers.get("gene", [""])[0]) for f in rec.features]
    assert order == [("source", ""), ("gene", "ND1"), ("CDS", "ND1"), ("gene", "trnY"), ("tRNA", "trnY")], order


def test_rrna_trimmed_for_inside_trna_near_the_start_up_to_250_nt():
    """Dev regression after review fix 7 (NC_023929.1 / NC_023931.1): MitoZ's l-rRNA began 199 nt before trnV, which
    therefore counted as 'inside' (loss 199 nt > 10 % of 1370). A tRNA-sized overhang at either end is a normal
    core error and is trimmed: the discarded side may be up to max(10 %, 250 nt)."""
    rec = _rrna_rec([("rRNA", 1001, 2370, "l-rRNA"), ("tRNA", 1011, 1200, "trnV(uac)")])
    rows = _curate_rec(rec)
    r = _by(rec, "rRNA", "rrnL")
    assert (int(r.location.start) + 1, int(r.location.end)) == (1201, 2370)
    assert any(x["rule"] == "R1" for x in rows) and not any(x["rule"] == "R3" for x in rows)


def test_rrna_overhang_is_trimmed_before_p1_pushes_a_cds_start():
    """Dev regression after review fix 12 (NC_060482.1 nad6 exact -> +9 nt, NC_045069.1 nad1 exact -> +30 nt): with the
    CDS loop before R1, the raw rRNA overhang still overlapped the CDS start and P1 moved the start downstream.
    Order is rRNA (R1) -> CDS -> rRNA again (so P3/O1 extensions are trimmed as well)."""
    from Bio.SeqFeature import FeatureLocation
    seq = ("CCCGGG" * 167)[:1000] + ("CCCGGG" * 84)[:500] + "ATG" + "GCC" * 32 + "ATG" + "GCC" * 60 + "TAA" + ("CCCGGG" * 300)[:1215]
    assert len(seq) == 3000 and seq[1500:1503] == "ATG" and seq[1599:1602] == "ATG"   # in-frame ATG at 1501 and 1600
    rec = _mkrec(seq, [("rRNA", FeatureLocation(1000, 1560, 1), "l-rRNA", {"product": ["16S ribosomal RNA"]}),
                       ("CDS", FeatureLocation(1500, 1785, 1), "ND6")])
    rows = _curate_rec(rec)
    cds, rrna = _by(rec, "CDS", "ND6"), _by(rec, "rRNA", "rrnL")
    assert int(cds.location.start) + 1 == 1501 and not any(r["rule"] == "P1" for r in rows), rows
    assert int(rrna.location.end) == 1500 and any(r["rule"] == "R1" for r in rows)


def test_minus_strand_frame_anchored_at_the_genome_end():
    """Dev regression NC_015829.1 (tardigrade) after review fix 6: nad4L/nad4 on the minus strand overlap by 7 nt in
    DIFFERENT frames (NCBI keeps the overlap; the frame of a minus-strand CDS is set by its genome END), but the
    frame test used the genome starts, called them 'same frame', and P1 moved the nad4 start 18 nt downstream.
    No P1 and no P2 for such a pair; a true same-frame overlap (multiple of 3) is still reported."""
    from Bio.SeqFeature import FeatureLocation
    from Bio.Seq import Seq
    nad4 = "ATG" + "GCC" * 5 + "ATG" + "GCC" * 94 + "T"      # 304 nt (T-tail stop, as NC_015829.1 nad4: length % 3 == 1);
    #                                                          in-frame ATG 18 nt downstream = where P1 would move the start
    nad4l = "ATG" + "GCC" * 60 + "TAA"                       # 186 nt
    # genome (minus strand genes): ... [nad4l rc] overlapping the last 7 nt of [nad4 rc]
    tx = nad4l[:-7] + nad4                                   # nad4l stop TAA + 4 nt overlap nad4 start -> 7 nt overlap
    rc = str(Seq(tx).reverse_complement())
    seq = ("CCCGGG" * 200)[:1000] + rc + ("CCCGGG" * 200)[:1000]
    L = len(seq)
    g_nad4_s, g_nad4_e = 1000, 1000 + len(nad4)              # nad4 occupies the 3' part of tx -> low genome coords
    g_nad4l_s, g_nad4l_e = 1000 + len(nad4) - 7, 1000 + len(tx)
    assert g_nad4l_s < g_nad4_e and (g_nad4l_s - g_nad4_s) % 3 == 0        # the OLD (start-based) test says 'same frame'
    assert (g_nad4l_e - g_nad4_e) % 3 != 0                                  # the end-based frame differs
    rec = _mkrec(seq, [("CDS", FeatureLocation(g_nad4_s, g_nad4_e, -1), "ND4"),
                       ("CDS", FeatureLocation(g_nad4l_s, g_nad4l_e, -1), "ND4L")])
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "ND4").location.end) == g_nad4_e, rows
    assert not any(r["rule"] in ("P1", "P2") for r in rows), rows
    assert L > 0
