"""batch_mitoz.resolve_tools / optional post-steps / final-step helpers (review B findings 6, 7, 10, 11, 14, 15, C8).

Run: python3 -m pytest -q tests/test_resolve_tools.py
"""
import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import batch_mitoz as bm  # noqa: E402
import mito_env  # noqa: E402


@pytest.fixture(autouse=True)
def restore_environ():
    """resolve_tools sets YZZ_* variables for the child processes."""
    old = os.environ.copy()
    yield
    os.environ.clear()
    os.environ.update(old)


def parse(tmp_path, *extra, inp=None):
    return bm.build_parser().parse_args(["run", "-i", str(inp or tmp_path), "-o", str(tmp_path / "out"), *extra])


# ---- B6: phylum from an explicit clade -------------------------------------------------------------------
@pytest.mark.parametrize("clade,phylum", [
    ("Mollusca", "Mollusca"), ("Annelida-segmented-worms", "Annelida"), ("Chordata", "Chordata"),
    ("Nemertea-ribbon-worms", "Nemertea"), ("Porifera-sponges", "Porifera"), ("Cnidaria", "Cnidaria"),
    ("Platyhelminthes-flatworms", "Platyhelminthes"), ("Chaetognatha-arrow-worms", "Chaetognatha")])
def test_explicit_clade_gives_the_phylum(tmp_path, clade, phylum):
    a = parse(tmp_path, "--clade", clade, "--genetic_code", "5")
    assert bm.resolve_tools(a) == "" and a.phylum == phylum


def test_animal_clade_leaves_phylum_unknown_and_explicit_phylum_wins(tmp_path):
    a = parse(tmp_path, "--clade", "Animal", "--genetic_code", "5")
    assert bm.resolve_tools(a) == "" and a.phylum is None
    a = parse(tmp_path, "--clade", "Animal", "--genetic_code", "5", "--phylum", "Ctenophora")
    assert bm.resolve_tools(a) == "" and a.phylum == "Ctenophora"
    a = parse(tmp_path, "--clade", "Mollusca", "--phylum", "Kinorhyncha")
    assert bm.resolve_tools(a) == "" and a.phylum == "Kinorhyncha"


def test_auto_clade_gives_the_phylum(tmp_path):
    r = SeqRecord(Seq("ATG" * 100), id="X.1", name="X", description="t")
    r.annotations.update(molecule_type="DNA", topology="circular", organism="G s",
                         taxonomy=["Eukaryota", "Metazoa", "Mollusca"])
    f = tmp_path / "m.gb"
    SeqIO.write(r, str(f), "genbank")
    a = parse(tmp_path, inp=f)
    assert bm.resolve_tools(a) == "" and (a.clade, a.phylum) == ("Mollusca", "Mollusca")


# ---- B10: the decisions are shown on the console --------------------------------------------------------
def test_clade_decision_is_printed(tmp_path, capsys):
    r = SeqRecord(Seq("ATG" * 100), id="X.1", name="X", description="t")
    r.annotations.update(molecule_type="DNA", topology="circular", organism="G s",
                         taxonomy=["Eukaryota", "Metazoa", "Mollusca"])
    f = tmp_path / "m.gb"
    SeqIO.write(r, str(f), "genbank")
    bm.resolve_tools(parse(tmp_path, inp=f))
    assert "[CLADE] auto: Mollusca -> --clade Mollusca --genetic_code 5" in capsys.readouterr().err


@pytest.mark.skipif(not (mito_env.mitoz_env() or mito_env.find("mitoz")), reason="MitoZ not installed here")
def test_env_decision_is_printed(tmp_path, capsys):
    bm.resolve_tools(parse(tmp_path, "--clade", "Mollusca"))
    assert "[ENV] MitoZ found automatically" in capsys.readouterr().err


# ---- B14: --core mitoz ------------------------------------------------------------------------------------
def test_core_mitoz_override_is_logged_and_documented(tmp_path, capsys):
    a = parse(tmp_path, "--clade", "Mollusca", "--core", "mitoz", "--rrna_hmm", "yes")
    assert bm.resolve_tools(a) == ""
    err = capsys.readouterr().err
    assert "[CORE]" in err and "--rrna_hmm yes" in err and "--rrna_hmm no" in err
    assert (a.trna, a.rescue, a.rrna_hmm) == ("mitoz", "no", "no")
    h = subprocess.run([sys.executable, str(HERE.parent / "batch_mitoz.py"), "run", "-h"], capture_output=True, text=True)
    assert "--rescue no --rrna_hmm no" in " ".join(h.stdout.split())


# ---- B7: optional steps never abort the run ---------------------------------------------------------------
def ns(**kw):
    base = dict(rescue="yes", rrna_hmm="yes", mitos_bin="/x", mitos_refdir="/y", genetic_code=5, phylum=None,
                rrna_hmm_lib=None)
    base.update(kw)
    return argparse.Namespace(**base)


def test_rrna_hmm_failure_is_logged_and_the_run_continues(tmp_path, monkeypatch, caplog):
    import mito_rrna_hmm

    def boom(*a, **k):
        raise subprocess.CalledProcessError(1, "nhmmer")
    monkeypatch.setattr(mito_rrna_hmm, "read_manifest", lambda lib: {"x": 1})
    monkeypatch.setattr(mito_env, "find", lambda *a, **k: "/bin/true")
    monkeypatch.setattr(mito_rrna_hmm, "search_dir", boom)
    with caplog.at_level(logging.ERROR):
        assert bm.run_rrna_hmm(ns(), tmp_path, tmp_path) is None
    assert any("[RRNA_HMM] failed" in m for m in caplog.messages)


def test_rescue_failure_is_logged_and_the_run_continues(tmp_path, monkeypatch, caplog):
    import mito_rescue

    def boom(*a, **k):
        raise RuntimeError("MITOS2 exploded")
    monkeypatch.setattr(mito_rescue, "rescue_dir", boom)
    with caplog.at_level(logging.ERROR):
        assert bm.run_rescue(ns(), tmp_path, tmp_path) is None
    assert any("[RESCUE] failed" in m and "MITOS2 exploded" in m for m in caplog.messages)


def test_wrapped_trna_failure_is_logged_and_the_run_continues(tmp_path, monkeypatch, caplog):
    import mito_trna

    def boom(*a, **k):
        raise ValueError("bad row")
    r = SeqRecord(Seq("ACGT" * 50), id="S1", description="")
    r.annotations["molecule_type"] = "DNA"
    SeqIO.write(r, str(tmp_path / "S1.gbf"), "genbank")
    (tmp_path / "wrap.tsv").write_text("S1\tY\tgta\tACGT\n")
    monkeypatch.setattr(mito_trna, "add_wrapped", boom)
    with caplog.at_level(logging.ERROR):
        assert bm.add_wrapped_trnas(tmp_path, tmp_path / "wrap.tsv") is None
    assert any("[TRNA_WRAP] failed" in m for m in caplog.messages)


# ---- B6/A: tell the user when the rRNA HMM step skipped records because no phylum was given --------------
@pytest.mark.parametrize("status", ["skipped: no library for this phylum", "skipped: no model for Kinorhyncha"])
def test_rrna_hmm_warns_when_nothing_could_be_searched(tmp_path, monkeypatch, caplog, capsys, status):
    import mito_rrna_hmm
    rows = [dict(file="S1.gbf", record="S1", phylum="", gene="s-rRNA", location="", status=status)]
    monkeypatch.setattr(mito_rrna_hmm, "read_manifest", lambda lib: {"x": 1})
    monkeypatch.setattr(mito_env, "find", lambda *a, **k: "/bin/true")
    monkeypatch.setattr(mito_rrna_hmm, "search_dir", lambda *a, **k: rows)
    with caplog.at_level(logging.WARNING):
        bm.run_rrna_hmm(ns(), tmp_path, tmp_path)
    assert any("--phylum" in m and "[RRNA_HMM]" in m for m in caplog.messages)
    assert "--phylum" in capsys.readouterr().err
    assert (tmp_path / "rrna_hmm_report.tsv").read_text().splitlines()[1].split("\t")[:6][-1] == status


def test_rrna_hmm_no_warning_when_models_were_used(tmp_path, monkeypatch, caplog, capsys):
    import mito_rrna_hmm
    rows = [dict(file="S1.gbf", record="S1", phylum="Mollusca", gene="s-rRNA", location="1..2", status="added")]
    monkeypatch.setattr(mito_rrna_hmm, "read_manifest", lambda lib: {"x": 1})
    monkeypatch.setattr(mito_env, "find", lambda *a, **k: "/bin/true")
    monkeypatch.setattr(mito_rrna_hmm, "search_dir", lambda *a, **k: rows)
    with caplog.at_level(logging.WARNING):
        bm.run_rrna_hmm(ns(), tmp_path, tmp_path)
    assert not any("--phylum" in m for m in caplog.messages) and "--phylum" not in capsys.readouterr().err


# ---- coordinator: report columns ---------------------------------------------------------------------------
def test_report_tables_carry_hint_and_coverage_columns(tmp_path, monkeypatch):
    import mito_rescue
    import mito_rrna_hmm
    monkeypatch.setattr(mito_rescue, "rescue_dir", lambda *a, **k: [
        dict(file="a.gbf", record="S1", gene="ND2", location="1..9", status="added", hint="check start")])
    bm.run_rescue(ns(), tmp_path, tmp_path)
    assert (tmp_path / "rescue_report.tsv").read_text().splitlines() == [
        "file\trecord\tgene\tlocation\tstatus\thint", "a.gbf\tS1\tND2\t1..9\tadded\tcheck start"]
    monkeypatch.setattr(mito_rescue, "rescue_dir", lambda *a, **k: [      # an older row without "hint"
        dict(file="a.gbf", record="S1", gene="ND2", location="1..9", status="added")])
    bm.run_rescue(ns(), tmp_path, tmp_path)
    assert (tmp_path / "rescue_report.tsv").read_text().splitlines()[1] == "a.gbf\tS1\tND2\t1..9\tadded\t"
    monkeypatch.setattr(mito_rrna_hmm, "read_manifest", lambda lib: {"x": 1})
    monkeypatch.setattr(mito_env, "find", lambda *a, **k: "/bin/true")
    monkeypatch.setattr(mito_rrna_hmm, "search_dir", lambda *a, **k: [
        dict(file="a.gbf", record="S1", phylum="Mollusca", gene="s-rRNA", location="1..9", status="added",
             coverage="0.93", note="check ends"),
        dict(file="a.gbf", record="S2", phylum="Mollusca", gene="l-rRNA", location="1..9", status="added (partial 5')")])
    bm.run_rrna_hmm(ns(), tmp_path, tmp_path)
    lines = (tmp_path / "rrna_hmm_report.tsv").read_text().splitlines()
    assert lines[0] == "file\trecord\tphylum\tgene\tlocation\tstatus\tcoverage\tnote"
    assert lines[1].endswith("\tadded\t0.93\tcheck ends") and lines[2].endswith("\tadded (partial 5')\t\t")


def test_added_count_includes_partial_statuses(tmp_path, monkeypatch, caplog):
    """Statuses like "added (partial 5')" are additions too."""
    import mito_rrna_hmm
    rows = [dict(file="a.gbf", record="S1", phylum="Mollusca", gene="s-rRNA", location="1..9", status="added (partial 5')"),
            dict(file="a.gbf", record="S2", phylum="Mollusca", gene="l-rRNA", location="1..9", status="rejected: short")]
    monkeypatch.setattr(mito_rrna_hmm, "read_manifest", lambda lib: {"x": 1})
    monkeypatch.setattr(mito_env, "find", lambda *a, **k: "/bin/true")
    monkeypatch.setattr(mito_rrna_hmm, "search_dir", lambda *a, **k: rows)
    with caplog.at_level(logging.INFO):
        bm.run_rrna_hmm(ns(), tmp_path, tmp_path)
    assert any("1 rRNA(s) added" in m and "S1:s-rRNA" in m and "S2" not in m for m in caplog.messages)


# ---- B11: nothing to re-annotate ----------------------------------------------------------------------------
def test_empty_phase2_list_starts_no_annotation_stage(tmp_path, monkeypatch):
    def boom(**k):
        raise AssertionError("run_annotation_stage must not be called for an empty list")
    monkeypatch.setattr(bm, "run_annotation_stage", boom)
    assert bm.run_phase2([], argparse.Namespace(), tmp_path, tmp_path) == {}


# ---- B15: body_source tag ------------------------------------------------------------------------------------
def test_body_source_tag_tells_auto_skip_from_phase2_failure():
    assert bm.fallback_src_tag("auto", "S1", {"S2"}) == "shifted_round1"          # auto decided: round 1 is fine
    assert bm.fallback_src_tag("auto", "S2", {"S2"}) == "phase1b_reoriented"       # re-annotation failed
    assert bm.fallback_src_tag("yes", "S1", {"S1"}) == "phase1b_reoriented"
    assert bm.fallback_src_tag("no", "S1", set()) == "phase1b_reoriented"


# ---- C8 / B9: tRNAs of samples finalised from fast round 1 --------------------------------------------------
def test_which_samples_get_the_sensitive_trna_pass(tmp_path):
    rows = [(tmp_path / "A.gbf", "phase2_mitoz"), (tmp_path / "B.gbf", "phase1_raw"),
            (tmp_path / "C.gbf", "phase1b_reoriented"), (tmp_path / "D.gbf", "phase1b_reoriented+partial_no_gene_merged"),
            (tmp_path / "E.gbf", "shifted_round1"), (tmp_path / "F.gbf", "phase2_mitoz+partial_no_gene_merged")]
    got = bm.phase1_fast_finalised(rows, phase1_fast=True)
    assert [p.name for p in got] == ["B.gbf", "C.gbf", "D.gbf"]
    assert bm.phase1_fast_finalised(rows, phase1_fast=False) == []


def test_sensitive_trna_pass_skipped_silently_without_mitos2(tmp_path, monkeypatch):
    def boom(*a, **k):
        raise AssertionError("no subprocess without MITOS2")
    monkeypatch.setattr(bm.subprocess, "run", boom)
    (tmp_path / "B.gbf").write_text("x")
    a = ns(mitos_bin=None, mitos_refdir=None, trna="mitos2")
    assert bm.resensitise_trnas(a, [tmp_path / "B.gbf"], tmp_path) == 0
    a = ns(trna="mitoz")                                    # MiTFi tRNAs requested: nothing is replaced
    assert bm.resensitise_trnas(a, [tmp_path / "B.gbf"], tmp_path) == 0
