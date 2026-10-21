"""Tests for benchmark/mito_bench_run.py (run tools on bare records, lift results to NCBI coordinates). No external tools.

Run: python3 -m pytest -q tests/test_mito_bench_run.py
"""
import json
import random
import sys
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "benchmark"))
import mito_bench_run as mr  # noqa: E402

L = 3000
random.seed(5)
NCBI = Seq("".join(random.choice("ACGT") for _ in range(L)))


def feat(loc, name="cox1", kind="gene"):
    return SeqFeature(loc, type=kind, qualifiers={"gene": [name]})


def rec(seq, feats):
    r = SeqRecord(seq, id="X")
    r.annotations["molecule_type"] = "DNA"
    r.features = feats
    return r


def lifted_seq(pred, lifted, i=0):
    return str(lifted.features[i].extract(NCBI)), str(pred.features[i].extract(pred.seq))


def test_lift_identity():
    p = rec(NCBI, [feat(FeatureLocation(100, 400, 1))])
    out = mr.lift_to_ncbi(p, NCBI)
    assert (int(out.features[0].location.start), int(out.features[0].location.end)) == (100, 400)
    assert out.seq == NCBI


def test_lift_rotation_origin_span():
    rot = 1000
    pseq = NCBI[rot:] + NCBI[:rot]                     # pred position i = ncbi (i + rot) % L
    # ncbi 900..1100 -> pred join(L-100..L, 0..100)
    p = rec(pseq, [feat(CompoundLocation([FeatureLocation(L - 100, L, 1), FeatureLocation(0, 100, 1)])),
                   feat(FeatureLocation(L - 1500, L - 900, 1), "nad1")])      # ncbi 2500..3000+100 spans origin
    out = mr.lift_to_ncbi(p, NCBI)
    loc = out.features[0].location
    assert (int(loc.start), int(loc.end), len(loc.parts)) == (900, 1100, 1)
    a, b = lifted_seq(p, out)
    assert a == b
    assert len(out.features[1].location.parts) == 2
    a, b = lifted_seq(p, out, 1)
    assert a == b


def test_lift_reverse_complement():
    pseq = (NCBI[500:] + NCBI[:500]).reverse_complement()
    p = rec(pseq, [feat(FeatureLocation(200, 800, 1)), feat(FeatureLocation(1000, 1300, -1), "rrnL", "rRNA")])
    out = mr.lift_to_ncbi(p, NCBI)
    assert out.features[0].location.strand == -1 and out.features[1].location.strand == 1
    for i in (0, 1):
        a, b = lifted_seq(p, out, i)
        assert a == b


def test_lift_length_mismatch():
    with pytest.raises(mr.LiftError):
        mr.lift_to_ncbi(rec(NCBI[:-5], []), NCBI)


ROW = {"accession": "NC_1.1", "genetic_code": "5", "mitoz_clade": "Annelida-segmented-worms", "topology": "circular"}


def test_build_cmd_unsupported(tmp_path):
    assert mr.build_cmd("mitoz", dict(ROW, mitoz_clade=""), tmp_path / "bare", tmp_path / "out", 2) is None
    assert mr.build_cmd("ours", dict(ROW, mitoz_clade=""), tmp_path / "bare", tmp_path / "out", 2) is None


def test_failed_tool_writes_status(tmp_path):
    exe = tmp_path / "fail.sh"
    exe.write_text("#!/bin/sh\necho boom >&2\nexit 1\n")
    exe.chmod(0o755)
    st = mr.run_steps([[str(exe)]], tmp_path / "o", timeout=10)
    assert st["status"] == "failed" and "boom" in st["error"]


def test_timeout_status(tmp_path):
    exe = tmp_path / "slow.sh"
    exe.write_text("#!/bin/sh\nsleep 5\n")
    exe.chmod(0o755)
    st = mr.run_steps([[str(exe)]], tmp_path / "o", timeout=1)
    assert st["status"] == "timeout"


def test_record_status_json(tmp_path):
    mr.write_status(tmp_path, {"tool": "mitoz", "accession": "NC_1.1", "status": "unsupported", "seconds": 0.0,
                               "step_seconds": {}, "error": ""})
    assert json.loads((tmp_path / "status.json").read_text())["status"] == "unsupported"


BED = "\n".join(["X\t26\t825\trrnS\t0.0\t+", "X\t822\t886\ttrnM_0(cat)\t1e-9\t+", "X\t900\t1500\tcox1-a\t1.0\t+",
                 "X\t1700\t2100\tcox1-b\t1.0\t+", "X\t2200\t2300\tOH_1\t9\t-", "X\t2300\t2600\tnad6\t1.0\t-"])


def test_mitos_to_record(tmp_path):
    (tmp_path / "result.bed").write_text(BED + "\n")
    r = mr.mitos_to_record(tmp_path, NCBI, "X")
    got = {(f.type, f.qualifiers["gene"][0]) for f in r.features if f.type != "gene"}
    assert got == {("rRNA", "rrnS"), ("tRNA", "trnM"), ("CDS", "cox1"), ("CDS", "nad6")}
    trna = next(f for f in r.features if f.type == "tRNA")
    assert trna.qualifiers["anticodon"] == ["cat"]
    cox1 = next(f for f in r.features if f.type == "CDS" and f.qualifiers["gene"] == ["cox1"])
    assert len(cox1.location.parts) == 2
    nad6 = next(f for f in r.features if f.type == "CDS" and f.qualifiers["gene"] == ["nad6"])
    assert nad6.location.strand == -1 and int(nad6.location.start) == 2300


def test_mitos2_outdir_exists(tmp_path, monkeypatch):
    """runmitos refuses a missing -o directory: run_one must create it (2026-10-20 smoke run)."""
    exe = tmp_path / "runmitos"
    exe.write_text('#!/bin/sh\nwhile [ "$1" != "-o" ]; do shift; done\n[ -d "$2" ] || { echo "no such directory $2" >&2; exit 2; }\n'
                   'printf "X\\t0\\t100\\tcox1\\t1\\t+\\n" > "$2/result.bed"\n')
    exe.chmod(0o755)
    monkeypatch.setattr(mr.mito_env, "find", lambda name, hint=None: str(exe))
    monkeypatch.setattr(mr.mito_env, "mitos_refdir", lambda a=None, b=None: str(tmp_path))
    (tmp_path / "gb").mkdir()
    from Bio import SeqIO
    SeqIO.write(rec(NCBI, []), str(tmp_path / "gb" / "NC_1.1.gb"), "genbank")
    st = mr.run_one("mitos2", ROW, tmp_path, 1, 30)
    assert st["status"] == "ok", st["error"]


def test_timeout_kills_children(tmp_path):
    """A timed-out step must not leave its child processes running (2026-10-20: orphan tblastn after timeout)."""
    pidfile = tmp_path / "child.pid"
    exe = tmp_path / "parent.sh"
    exe.write_text(f"#!/bin/sh\nsleep 30 &\necho $! > {pidfile}\nwait\n")
    exe.chmod(0o755)
    st = mr.run_steps([[str(exe)]], tmp_path / "o", timeout=1)
    assert st["status"] == "timeout"
    import os
    import time
    pid = int(pidfile.read_text())
    time.sleep(0.5)
    try:
        os.kill(pid, 0)
        alive = Path(f"/proc/{pid}/stat").read_text().split()[2] != "Z"
    except ProcessLookupError:
        alive = False
    assert not alive


def test_mitos_copies_kept(tmp_path):
    """Review C1: MITOS2 names extra copies nad5_0 / nad5_1; both must stay nad5."""
    (tmp_path / "result.bed").write_text("X\t100\t700\tnad5_0\t1\t+\nX\t1000\t1600\tnad5_1\t1\t+\n")
    r = mr.mitos_to_record(tmp_path, NCBI, "X")
    assert [f.qualifiers["gene"][0] for f in r.features if f.type == "CDS"] == ["nad5", "nad5"]


def test_mitos_fragments_in_letter_order(tmp_path):
    """Review C2: origin-spanning gene as cox2-b 1..268 + cox2-a 2700..3000 -> join(a, b)."""
    (tmp_path / "result.bed").write_text("X\t0\t268\tcox2-b\t1\t+\nX\t2700\t3000\tcox2-a\t1\t+\n")
    r = mr.mitos_to_record(tmp_path, NCBI, "X")
    cds = next(f for f in r.features if f.type == "CDS")
    assert [int(p.start) for p in cds.location.parts] == [2700, 0]
    assert str(cds.extract(NCBI)) == str(NCBI[2700:] + NCBI[:268])


def test_mitos_start_after_end_wraps(tmp_path):
    """Review I5: MITOS circular features may be written with start > end."""
    (tmp_path / "result.bed").write_text("X\t2900\t100\tcox1\t1\t+\n")
    r = mr.mitos_to_record(tmp_path, NCBI, "X")
    cds = next(f for f in r.features if f.type == "CDS")
    assert str(cds.extract(NCBI)) == str(NCBI[2900:] + NCBI[:100])


def test_run_one_unexpected_error_is_failed(tmp_path, monkeypatch):
    """Review I3: an unexpected exception while collecting must give status failed, not abort the run."""
    exe = tmp_path / "runmitos"
    exe.write_text('#!/bin/sh\nexit 0\n')
    exe.chmod(0o755)
    monkeypatch.setattr(mr.mito_env, "find", lambda name, hint=None: str(exe))
    monkeypatch.setattr(mr.mito_env, "mitos_refdir", lambda a=None, b=None: str(tmp_path))
    monkeypatch.setattr(mr, "collect", lambda *a: (_ for _ in ()).throw(KeyError("boom")))
    (tmp_path / "gb").mkdir()
    from Bio import SeqIO
    SeqIO.write(rec(NCBI, []), str(tmp_path / "gb" / "NC_1.1.gb"), "genbank")
    st = mr.run_one("mitos2", ROW, tmp_path, 1, 30)
    assert st["status"] == "failed" and "boom" in st["error"]


def test_recollect_reuses_tool_output(tmp_path, monkeypatch):
    """Review I6: --recollect rebuilds final.gbf from kept tool output without running the tool again."""
    out = tmp_path / "runs" / "mitos2" / "NC_1.1"
    (out / "mitos").mkdir(parents=True)
    (out / "mitos" / "result.bed").write_text("X\t100\t700\tnad5_0\t1\t+\n")
    mr.write_status(out, {"tool": "mitos2", "accession": "NC_1.1", "status": "ok", "seconds": 5.0,
                          "step_seconds": {}, "error": ""})
    (tmp_path / "gb").mkdir()
    from Bio import SeqIO
    SeqIO.write(rec(NCBI, []), str(tmp_path / "gb" / "NC_1.1.gb"), "genbank")
    monkeypatch.setattr(mr, "run_steps", lambda *a, **k: pytest.fail("tool must not run"))
    st = mr.run_one("mitos2", ROW, tmp_path, 1, 30, recollect=True)
    assert st["status"] == "ok" and st["seconds"] == 5.0
    genes = [f.qualifiers["gene"][0] for f in SeqIO.read(str(out / "final.gbf"), "genbank").features if f.type == "CDS"]
    assert genes == ["nad5"]


def test_recollect_retries_collect_failures(tmp_path, monkeypatch):
    """2026-10-20 dev run: 10 MITOS2 records failed in collect/lift under the old parser (start > end).
    --recollect must rebuild them too, not only 'ok' runs."""
    out = tmp_path / "runs" / "mitos2" / "NC_1.1"
    (out / "mitos").mkdir(parents=True)
    (out / "mitos" / "result.bed").write_text("X\t2900\t100\tcox1\t1\t+\n")
    mr.write_status(out, {"tool": "mitos2", "accession": "NC_1.1", "status": "failed", "seconds": 5.0,
                          "step_seconds": {}, "error": "collect/lift: End location (100) must be greater"})
    (tmp_path / "gb").mkdir()
    from Bio import SeqIO
    SeqIO.write(rec(NCBI, []), str(tmp_path / "gb" / "NC_1.1.gb"), "genbank")
    monkeypatch.setattr(mr, "run_steps", lambda *a, **k: pytest.fail("tool must not run"))
    st = mr.run_one("mitos2", ROW, tmp_path, 1, 30, recollect=True)
    assert st["status"] == "ok" and st["error"] == "" and st["seconds"] == 5.0


def test_recollect_keeps_tool_failures(tmp_path):
    """A tool that itself failed has no output to recollect: status stays failed."""
    out = tmp_path / "runs" / "mitos2" / "NC_1.1"
    out.mkdir(parents=True)
    mr.write_status(out, {"tool": "mitos2", "accession": "NC_1.1", "status": "failed", "seconds": 5.0,
                          "step_seconds": {}, "error": "runmitos: boom"})
    assert mr.run_one("mitos2", ROW, tmp_path, 1, 30, recollect=True)["error"] == "runmitos: boom"


def test_mitom_runs_every_phylum(tmp_path):
    """MitoM (hybrid core) covers phyla MitoZ lacks via the extended clades: never 'unsupported'."""
    row = dict(ROW, mitoz_clade="", phylum="Cnidaria", genetic_code="4")
    cmd = mr.build_cmd("mitom", row, tmp_path / "bare", tmp_path / "out", 2)
    assert cmd[0][cmd[0].index("--clade") + 1] == "Cnidaria" and cmd[0][cmd[0].index("--genetic_code") + 1] == "4"
    cmd = mr.build_cmd("mitom", dict(row, phylum="Rotifera"), tmp_path / "bare", tmp_path / "out", 2)
    assert cmd[0][cmd[0].index("--clade") + 1] == "Animal"
    cmd = mr.build_cmd("mitom", ROW, tmp_path / "bare", tmp_path / "out", 2)
    assert cmd[0][cmd[0].index("--clade") + 1] == "Annelida-segmented-worms"


def test_mitom_auto_differs_only_in_reannotation(tmp_path):
    """A/B of --if_reannotation_2 auto: same command as mitom plus the option; output collected the same way."""
    a = mr.build_cmd("mitom", ROW, tmp_path / "bare", tmp_path / "out", 2)[0]
    b = mr.build_cmd("mitom_auto", ROW, tmp_path / "bare", tmp_path / "out", 2)[0]
    assert b == a + ["--if_reannotation_2", "auto"]
    assert "mitom_auto" in mr.TOOLS
