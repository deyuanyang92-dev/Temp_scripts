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
