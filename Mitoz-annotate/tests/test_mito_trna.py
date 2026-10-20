"""Tests for mito_trna.py (multi-tool tRNA annotation, GeSeq-style review + selected GenBank).

Unit tests need no external tools; the integration test runs MITOS2/ARWEN/ARAGORN/tRNAscan-SE on
Orbicella annularis NC_007224 (bare sequence, NCBI annotates only trnM and trnW) and is skipped when
the tools are not installed. Run: python3 -m pytest -q tests/test_mito_trna.py
"""
import csv
import os
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_trna as mt  # noqa: E402

TRNA_BIN = Path(os.environ.get("MITO_TRNA_BIN", "/root/miniconda3/envs/mitotrna/bin"))
MITOS_BIN = Path(os.environ.get("MITOS_BIN", "/root/miniconda3/envs/mitos/bin"))
MITOS_REF = os.environ.get("MITOS_REFDIR", "")


def test_laslett_origin_spanning_unwrapped():
    text = "1   mtRNA-Met                c[16138,73]     35      (cat)\n2   mtRNA-Trp   [13399,13469]  34  (tca)\n"
    p = mt.unwrap(mt.parse_laslett(text, "ARWEN"), 16138)
    assert (p[0].aa, p[0].start, p[0].end, p[0].strand, p[0].anticodon) == ("M", 16138, 16211, -1, "cat")
    assert (p[1].aa, p[1].start, p[1].end, p[1].strand) == ("W", 13399, 13469, 1)


def test_mitos_bed_keeps_numbered_copies():
    # 2026-10-20: trnM_0(cat) / trnY_1(gta) (extra copies, Lottia digitalis NC_007782) were silently dropped
    bed = ("x\t5774\t5841\ttrnM_0(cat)\t7.2e-08\t+\nx\t6857\t6927\ttrnY_1(gta)\t8.4e-10\t-\n"
           "x\t100\t165\ttrnL1(tag)\t1e-9\t+\nx\t15508\t15575\ttrnQ(---)\t4.8e-4\t+\nx\t200\t260\tcox1\t0\t+\n")
    p = mt.parse_mitos_bed(bed)
    assert [(x.aa, x.anticodon, x.start, x.end, x.strand) for x in p] == [
        ("M", "cat", 5775, 5841, 1), ("Y", "gta", 6858, 6927, -1), ("L", "tag", 101, 165, 1),
        ("Q", "", 15509, 15575, 1)]   # trnQ(---): anticodon unreadable, tRNA kept


def test_origin_spanning_groups_with_unwrapped_copy():
    """MITOS2 1..72 and ARWEN 16138..16211 (= 16138,1..73) are the same tRNA on a circular genome."""
    mt.CIRC_L = 16138
    try:
        g = mt.group([mt.Pred("MITOS2", 1, 72, 1, "M", "cat"), mt.Pred("ARWEN", 16138, 16211, 1, "M", "cat")])
    finally:
        mt.CIRC_L = 0
    assert len(g) == 1 and mt.families(g[0]) == ["ARWEN/ARAGORN (heuristic)", "MITOS2 (covariance models, metazoan mt tRNAs)"]


def test_arwen_aragorn_count_as_one_family():
    g = mt.group([mt.Pred("ARWEN", 100, 170, 1, "F", "gaa"), mt.Pred("ARAGORN", 100, 170, 1, "F", "gaa")])
    assert len(g) == 1 and len(mt.families(g[0])) == 1


def test_different_anticodon_not_grouped():
    g = mt.group([mt.Pred("MITOS2", 100, 165, 1, "L", "tag"), mt.Pred("ARWEN", 100, 165, 1, "L", "taa")])
    assert len(g) == 2 and mt.conflicts(g) == {1: [2], 2: [1]}


def test_leu_ser_names_by_anticodon():
    assert [mt.trna_name("L", "tag"), mt.trna_name("L", "taa"), mt.trna_name("S", "tct"),
            mt.trna_name("S", "gct"), mt.trna_name("S", "tga"), mt.trna_name("M", "cat")] == \
        ["trnL1", "trnL2", "trnS1", "trnS1", "trnS2", "trnM"]


def test_feature_across_origin_is_join():
    f = mt.feature(mt.Pred("ARWEN", 16138, 16211, 1, "M", "cat"), "cat", "x", 16138)
    assert str(f.location) == "join{[16137:16138](+), [0:73](+)}"
    assert f.qualifiers["product"] == ["tRNA-Met"]


@pytest.mark.skipif(not ((TRNA_BIN / "arwen").exists() and (MITOS_BIN / "runmitos").exists() and MITOS_REF),
                    reason="tRNA tools / MITOS2 reference (MITOS_REFDIR) not available")
def test_orbicella_selected_matches_ncbi(tmp_path):
    """Anthozoan with 2 tRNAs: ARWEN/ARAGORN report several extra tRNAs (kept in the review GenBank
    only); the selected GenBank keeps exactly trnM and trnW, as NCBI NC_007224."""
    subprocess.run([sys.executable, "-W", "ignore", str(HERE.parent / "mito_trna.py"), "-i", str(HERE / "data" / "NC_007224_bare.gbf"),
                    "-o", str(tmp_path), "--genetic_code", "4", "--tools", "mitos2,arwen,aragorn,trnascan",
                    "--bin", str(TRNA_BIN), "--mitos_bin", str(MITOS_BIN), "--mitos_refdir", MITOS_REF],
                   check=True, capture_output=True)
    sel = SeqIO.read(str(tmp_path / "selected" / "NC_007224_bare.gbf"), "genbank")
    assert sorted(f.qualifiers["product"][0] for f in sel.features if f.type == "tRNA") == ["tRNA-Met", "tRNA-Trp"]
    rev = SeqIO.read(str(tmp_path / "review" / "NC_007224_bare.gbf"), "genbank")
    assert sum(f.type == "tRNA" for f in rev.features) > 2
    rows = list(csv.DictReader((tmp_path / "trna_report.tsv").open(), delimiter="\t"))
    met = next(r for r in rows if r["tRNA"] == "trnM")
    assert (met["start"], met["end"], met["boundary_source"]) == ("1", "72", "MITOS2")


def test_mitos2_failure_keeps_input_trnas(tmp_path):
    """--select mitos2 must never delete the MitoZ tRNAs because MITOS2 crashed (no result.bed)."""
    fake = tmp_path / "bin"
    fake.mkdir()
    for n in ("runmitos", "python"):
        (fake / n).write_text("#!/bin/sh\nexit 1\n")
        (fake / n).chmod(0o755)
    subprocess.run([sys.executable, "-W", "ignore", str(HERE.parent / "mito_trna.py"), "-i", str(HERE / "data" / "KY753833_mitoz.gbf"),
                    "-o", str(tmp_path / "out"), "--select", "mitos2", "--mitos_bin", str(fake), "--mitos_refdir", str(tmp_path)],
                   check=True, capture_output=True)
    before = SeqIO.read(str(HERE / "data" / "KY753833_mitoz.gbf"), "genbank")
    after = SeqIO.read(str(tmp_path / "out" / "selected" / "KY753833_mitoz.gbf"), "genbank")
    assert [str(f.location) for f in after.features] == [str(f.location) for f in before.features]


def test_select_mitos2_without_mitos2_refuses(tmp_path):
    r = subprocess.run([sys.executable, "-W", "ignore", str(HERE.parent / "mito_trna.py"), "-i", str(HERE / "data" / "KY753833_mitoz.gbf"),
                        "-o", str(tmp_path / "out"), "--select", "mitos2", "--mitos_bin", str(tmp_path)], capture_output=True)
    assert r.returncode == 2 and not (tmp_path / "out").exists()


@pytest.mark.skipif(not ((MITOS_BIN / "runmitos").exists() and MITOS_REF), reason="MITOS2 / MITOS_REFDIR not available")
def test_mitos2_replaces_mitoz_trnas(tmp_path):
    """MitoZ annotation of Aphrodita KY753833: tRNAs (and their gene features) replaced by MITOS2;
    all PCG/rRNA features untouched."""
    src = HERE / "data" / "KY753833_mitoz.gbf"
    subprocess.run([sys.executable, "-W", "ignore", str(HERE.parent / "mito_trna.py"), "-i", str(src), "-o", str(tmp_path),
                    "--select", "mitos2", "--genetic_code", "5", "--mitos_bin", str(MITOS_BIN), "--mitos_refdir", MITOS_REF],
                   check=True, capture_output=True)
    before = SeqIO.read(str(src), "genbank")
    after = SeqIO.read(str(tmp_path / "selected" / "KY753833_mitoz.gbf"), "genbank")
    trna = [f for f in after.features if f.type == "tRNA"]
    assert len(trna) == 22 and all("boundaries from MITOS2" in f.qualifiers["note"][0] for f in trna)
    genes = [f.qualifiers["gene"][0] for f in after.features if f.type == "gene"]
    assert len(genes) == len(set(genes)) == sum(f.type == "gene" for f in before.features)
    keep = lambda r: sorted(str(f.location) for f in r.features if f.type in ("CDS", "rRNA"))  # noqa: E731
    assert keep(after) == keep(before)


def test_select_mitos2_refuses_unoriented_input(tmp_path):
    """Replacement must run on the final reoriented GenBank: FASTA or round-1 (01.mitoz_anno1) refused."""
    d = tmp_path / "01.mitoz_anno1"
    d.mkdir()
    (d / "x.gbf").write_bytes((HERE / "data" / "KY753833_mitoz.gbf").read_bytes())
    (tmp_path / "x.fasta").write_text(">x\nACGT\n")
    for src in (d, tmp_path / "x.fasta"):
        r = subprocess.run([sys.executable, "-W", "ignore", str(HERE.parent / "mito_trna.py"), "-i", str(src), "-o",
                            str(tmp_path / "out"), "--select", "mitos2", "--mitos_bin", str(MITOS_BIN),
                            "--mitos_refdir", str(tmp_path)], capture_output=True, text=True)
        if not (MITOS_BIN / "runmitos").exists():
            pytest.skip("MITOS2 not available")
        assert r.returncode == 2 and "reoriented" in r.stderr
