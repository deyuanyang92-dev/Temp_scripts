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
    sel = SeqIO.read(str(tmp_path / "NC_007224_bare.selected.gbf"), "genbank")
    assert sorted(f.qualifiers["product"][0] for f in sel.features if f.type == "tRNA") == ["tRNA-Met", "tRNA-Trp"]
    rev = SeqIO.read(str(tmp_path / "NC_007224_bare.review.gbf"), "genbank")
    assert sum(f.type == "tRNA" for f in rev.features) > 2
    rows = list(csv.DictReader((tmp_path / "trna_report.tsv").open(), delimiter="\t"))
    met = next(r for r in rows if r["tRNA"] == "trnM")
    assert (met["start"], met["end"], met["boundary_source"]) == ("1", "72", "MITOS2")
