"""Tests for benchmark/mito_bench_eval.py (agreement with the NCBI annotation, obvious errors). No external tools.

Run: python3 -m pytest -q tests/test_mito_bench_eval.py
"""
import sys
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "benchmark"))
import mito_bench_eval as me  # noqa: E402


def f(kind, s, e, strand=1, **q):
    return SeqFeature(FeatureLocation(s, e, strand), type=kind, qualifiers={k: [v] for k, v in q.items()})


def rec(feats, seq="A" * 5000):
    r = SeqRecord(Seq(seq), id="X")
    r.features = feats
    return r


def classes(truth_feats, pred_feats):
    rows = me.match(me.norm_features(rec(truth_feats)), me.norm_features(rec(pred_feats)))
    return sorted((r["name"], r["class"]) for r in rows)


def test_norm_names():
    fs = me.norm_features(rec([f("CDS", 0, 99, gene="COX1a"), f("CDS", 200, 299, product="cytochrome b"),
                               f("rRNA", 400, 999, product="16S ribosomal RNA"),
                               f("tRNA", 1000, 1066, product="tRNA-Leu", anticodon="(pos:1030..1032,aa:Leu,seq:tag)")]))
    assert [(x.kind, x.name, x.anticodon) for x in fs] == [("PCG", "cox1", ""), ("PCG", "cob", ""),
                                                           ("rRNA", "rrnL", ""), ("tRNA", "L", "tag")]


def test_exact_and_le3():
    t = [f("CDS", 100, 1600, gene="cox1")]
    assert classes(t, [f("CDS", 100, 1600, gene="COX1")]) == [("cox1", "exact")]
    assert classes(t, [f("CDS", 102, 1600, gene="cox1")]) == [("cox1", "le3")]
    assert classes(t, [f("CDS", 110, 1600, gene="cox1")]) == [("cox1", "gt3")]


def test_strand_error():
    assert classes([f("CDS", 100, 1600, gene="cox1")], [f("CDS", 100, 1600, -1, gene="cox1")]) == [("cox1", "strand_error")]


def test_missed_and_fp():
    got = classes([f("CDS", 100, 1600, gene="cox1")], [f("CDS", 3000, 3500, gene="nad6")])
    assert got == [("cox1", "missed"), ("nad6", "fp")]


def test_trnL_name_error():
    t = [f("tRNA", 1000, 1066, product="tRNA-Leu", anticodon="(pos:1030..1032,aa:Leu,seq:tag)")]
    p = [f("tRNA", 1000, 1066, gene="trnL2", product="tRNA-Leu", note="anticodon: taa")]
    assert classes(t, p) == [("L", "name_error")]


def test_no_anticodon_unverifiable():
    t = [f("tRNA", 1000, 1066, product="tRNA-Ser")]
    p = [f("tRNA", 1000, 1066, gene="trnS2(tga)")]
    assert classes(t, p) == [("S", "name_unverifiable")]


CODE5 = 5


def cds_rec(seq, **q):
    return rec([f("CDS", 0, len(seq), gene="cox3", **q)], seq=seq + "A" * 100)


def test_internal_stop():
    checks = me.pcg_checks(cds_rec("ATG" + "TAA" + "GCC" * 10 + "TAA"), CODE5)
    assert [c["check"] for c in checks] == ["internal_stop"]


def test_incomplete_stop_ok():
    assert me.pcg_checks(cds_rec("ATG" + "GCC" * 10 + "T"), CODE5) == []
    assert me.pcg_checks(cds_rec("ATG" + "GCC" * 10 + "TA"), CODE5) == []


def test_bad_start():
    assert [c["check"] for c in me.pcg_checks(cds_rec("CCC" + "GCC" * 10 + "TAA"), CODE5)] == ["bad_start"]


def test_unsupported_excluded():
    obv = me.obvious_errors([], [], {"status": "unsupported"})
    assert obv == []


def test_failed_is_obvious_error():
    obv = me.obvious_errors([], [], {"status": "failed", "error": "boom"})
    assert [o["error"] for o in obv] == ["failed"]


def test_obvious_from_rows():
    rows = [{"kind": "PCG", "name": "nad5", "class": "missed", "delta": ""},
            {"kind": "tRNA", "name": "L", "class": "name_error", "delta": ""},
            {"kind": "PCG", "name": "cox1", "class": "gt3", "delta": 45},
            {"kind": "PCG", "name": "cox2", "class": "gt3", "delta": 12},
            {"kind": "tRNA", "name": "W", "class": "missed", "delta": ""}]
    obv = me.obvious_errors(rows, [{"gene": "cox3", "check": "internal_stop"}], {"status": "ok"})
    assert sorted(o["error"] for o in obv) == ["boundary_gt30", "internal_stop", "missing_gene", "name_error"]


def test_duplicate_gene():
    p = me.norm_features(rec([f("CDS", 100, 600, gene="nad1"), f("CDS", 2000, 2500, gene="nad1")]))
    assert me.duplicates(p) == ["nad1"]
