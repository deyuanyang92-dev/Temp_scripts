"""Regression tests for blast2_metadatav2.12.py with a real local BLAST+ database.

Run: mamba run -n blast python -m pytest -q tests/
Needs blastn and makeblastdb on PATH (or BLAST_BIN=<dir>).
"""
import csv
import importlib.util
import io
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
SCRIPT = HERE.parent / "blast2_metadatav2.12.py"
BLAST_BIN = os.environ.get("BLAST_BIN") or str(Path(sys.executable).parent)
ENV = {**os.environ, "PATH": BLAST_BIN + os.pathsep + os.environ.get("PATH", "")}

pytestmark = pytest.mark.skipif(
    shutil.which("blastn", path=ENV["PATH"]) is None, reason="BLAST+ not installed"
)


def load_module():
    spec = importlib.util.spec_from_file_location("b2m", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["b2m"] = mod  # dataclasses resolves annotations via sys.modules
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def db(tmp_path_factory):
    d = tmp_path_factory.mktemp("db")
    subprocess.run(
        ["makeblastdb", "-in", str(DATA / "mito_db.fasta"), "-dbtype", "nucl",
         "-parse_seqids", "-out", str(d / "mito")],
        check=True, env=ENV, capture_output=True,
    )
    return d / "mito"


def run(args, cwd):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args], cwd=cwd, env=ENV,
        capture_output=True, text=True,
    )


def read_csv(p):
    with open(p, encoding="utf-8-sig", newline="") as f:  # results CSV has a BOM for Excel
        return list(csv.DictReader(f))


def local_args(query, out, db, *extra):
    return ["-i", str(query), "-o", str(out), "-bm", "local", "--blast_db_nt", str(db),
            "--no_ncbi_metadata", "--batch_size", "2", *extra]


def test_topn_counts_subjects_not_hsps(db, tmp_path):
    r = run(local_args(DATA / "query.fasta", tmp_path / "o", db,
                       "--no_smk", "--mode", "blast2topN", "--topN", "3"), tmp_path)
    assert r.returncode == 0, r.stdout[-2000:]
    rows = read_csv(tmp_path / "o" / "query_results.csv")
    per_q = {}
    for row in rows:
        per_q.setdefault(row["qseqid"], []).append(row["sseqid"])
    # LC670763 COX1 spans the origin -> 2 HSPs to the same subject in the raw table
    raw = (tmp_path / "o" / "blast_results" / "query_blast.tsv").read_text().splitlines()
    assert sum(1 for x in raw if x.startswith("LC670763.1_COX1\t")) >= 2
    for q, subs in per_q.items():
        assert len(subs) == len(set(subs)), (q, subs)
        assert q.split("_COX")[0] in subs[0], (q, subs)  # self hit ranks first


def test_bitscore_default_beats_short_perfect_hit():
    m = load_module()
    rows = [
        {"qseqid": "q", "sseqid": "short", "pident": 100.0, "qcovs": 5, "evalue": 1e-8, "bitscore": 56.0},
        {"qseqid": "q", "sseqid": "full", "pident": 99.0, "qcovs": 100, "evalue": 0.0, "bitscore": 1090.0},
    ]
    args = m.parse_arguments(["-i", "x", "--no_smk"])
    assert args.sort_by == "bitscore"
    assert m.get_topN_hits_py(rows, 1, args.sort_by)[0]["sseqid"] == "full"
    assert m.get_topN_hits_py(rows, 1, "identity")[0]["sseqid"] == "short"


def test_failed_batch_fails_file_and_resume_retries_only_it(db, tmp_path):
    fake = tmp_path / "flaky_blastn"
    fake.write_text(
        "#!/bin/sh\n"
        'case "$*" in *batch_00001__*) echo boom >&2; exit 1;; esac\n'
        f'exec {shutil.which("blastn", path=ENV["PATH"])} "$@"\n'
    )
    fake.chmod(0o755)
    out = tmp_path / "o"
    base = local_args(DATA / "query.fasta", out, db, "--no_smk", "--mode", "blast2topN",
                      "--retries", "1")
    r = run(base + ["--blast_exe", str(fake)], tmp_path)
    assert r.returncode != 0
    assert not (out / "query_completed.txt").exists()
    assert not (out / "query_results.csv").exists()
    assert not (out / "blast_results" / "query_blast.tsv").exists()
    batch_dir = out / "blast_results" / "query"
    good = batch_dir / "batch_00000_blast.tsv"
    assert good.stat().st_size > 0
    assert not (batch_dir / "batch_00001_blast.tsv").exists()
    assert not list(batch_dir.glob("*.part"))
    mtime = good.stat().st_mtime_ns

    r = run(base + ["--resume"], tmp_path)
    assert r.returncode == 0, r.stdout[-2000:]
    assert good.stat().st_mtime_ns == mtime  # finished batch not rerun
    qs = {row["qseqid"] for row in read_csv(out / "query_results.csv")}
    n_query = (DATA / "query.fasta").read_text().count(">")
    assert len(qs) == n_query


def test_resume_ignores_leftover_part_file(db, tmp_path):
    out = tmp_path / "o"
    batch_dir = out / "blast_results" / "query"
    batch_dir.mkdir(parents=True)
    (batch_dir / "batch_00000_blast.tsv.part").write_text("truncated\tline\n")
    r = run(local_args(DATA / "query.fasta", out, db, "--no_smk", "--mode", "blast2topN",
                       "--resume"), tmp_path)
    assert r.returncode == 0, r.stdout[-2000:]
    assert "truncated" not in (out / "query_results.csv").read_text()


def test_qblast_columns_match_blastplus(db, tmp_path):
    """Feed the qblast parser BLAST XML from the same local search and compare to outfmt 6."""
    m = load_module()
    m.BLAST_LIMITER = m._RateLimiter(1000, 1.0)
    query = DATA / "query.fasta"
    xml = subprocess.run(
        ["blastn", "-query", str(query), "-db", str(db), "-outfmt", "5", "-evalue", "1e-5"],
        check=True, env=ENV, capture_output=True, text=True,
    ).stdout
    tab = subprocess.run(
        ["blastn", "-query", str(query), "-db", str(db), "-evalue", "1e-5",
         "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"],
        check=True, env=ENV, capture_output=True, text=True,
    ).stdout
    # outfmt 5 has one <BlastOutput> per run; split it per query for the per-sequence qblast calls
    from Bio.Blast import NCBIXML
    per_query = {rec.query.split()[0]: rec for rec in NCBIXML.parse(io.StringIO(xml))}

    def fake_qblast(program, database, sequence, **kw):
        for qid, seq in seqs.items():
            if seq == sequence:
                return io.StringIO(qid)
        raise AssertionError("unknown query")

    from Bio import SeqIO
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(str(query), "fasta")}
    m.NCBIWWW.qblast = fake_qblast
    m.NCBIXML.parse = lambda h: [per_query[h.read()]]

    args = m.parse_arguments(["-i", str(query), "--no_smk", "-bm", "online",
                              "--online_engine", "biopython", "--online_delay_sec", "0"])
    batch = m.BatchInfo("batch_00000", query)
    _, ok, out_file, cols = m.run_biopython_qblast_one_batch(batch, "nucleotide", args, tmp_path, m.BLAST_HEADERS)
    assert ok
    got = m.parse_blast_results_py(out_file, cols)

    def key(r):
        return (r["qseqid"], m.accession_base(m.parse_accession_from_sseqid(str(r["sseqid"]))),
                int(r["qstart"]), int(r["qend"]), int(r["sstart"]), int(r["send"]))

    want = {}
    for line in tab.splitlines():
        p = line.split("\t")
        r = dict(zip(["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                      "qend", "sstart", "send", "evalue", "bitscore", "qcovs"], p))
        want[key(r)] = r
    assert len(got) == len(want)
    gaps_seen = False
    for r in got:
        w = want[key(r)]
        assert int(r["length"]) == int(w["length"])
        assert int(r["mismatch"]) == int(w["mismatch"]), (r, w)
        assert int(r["gapopen"]) == int(w["gapopen"]), (r, w)
        assert abs(float(r["pident"]) - float(w["pident"])) < 0.01
        assert abs(float(r["bitscore"]) - float(w["bitscore"])) < 1.0, (r, w)
        assert int(r["qcovs"]) == int(w["qcovs"]), (r, w)
        gaps_seen |= int(w["gapopen"]) > 0
    assert gaps_seen  # indel query exercised gap counting


def test_topn_only_from_existing_blast_table(db, tmp_path):
    r = run(local_args(DATA / "query.fasta", tmp_path / "a", db, "--no_smk", "--mode", "blast-only"), tmp_path)
    assert r.returncode == 0, r.stdout[-2000:]
    table = tmp_path / "a" / "blast_results" / "query_blast.tsv"
    r = run(["-i", str(table), "-o", str(tmp_path / "b"), "--no_smk", "--mode", "topN-only", "--topN", "1"],
            tmp_path)
    assert r.returncode == 0, r.stdout[-2000:]
    rows = read_csv(tmp_path / "b" / "query_results.csv")
    assert len(rows) == len({r["qseqid"] for r in rows}) == (DATA / "query.fasta").read_text().count(">")


def test_snakemake_scheduler_end_to_end(db, tmp_path):
    if shutil.which("snakemake", path=ENV["PATH"]) is None:
        pytest.skip("snakemake not installed")
    q = tmp_path / "in"
    q.mkdir()
    shutil.copy(DATA / "query.fasta", q / "s1.fasta")
    shutil.copy(DATA / "query.fasta", q / "s2.fasta")
    rel_db = os.path.relpath(db, tmp_path)  # child jobs run in another cwd; relative DB must still work
    r = run(local_args("in", "o", rel_db, "--smk", "--mode", "blast2topN", "--threads", "2"), tmp_path)
    assert r.returncode == 0, (r.stdout + r.stderr)[-3000:]
    for s in ("s1", "s2"):
        assert (tmp_path / "o" / s / f"{s}_completed.txt").exists()
        assert read_csv(tmp_path / "o" / s / f"{s}_results.csv")


# ── headers ──────────────────────────────────────────────────────────────

def blast_table(db, tmp_path, outfmt, name):
    out = tmp_path / name
    subprocess.run(["blastn", "-query", str(DATA / "query.fasta"), "-db", str(db), "-outfmt", outfmt,
                    "-evalue", "1e-5", "-out", str(out)], check=True, env=ENV, capture_output=True)
    return out


def add_header(table, tmp_path, *extra):
    r = run(["-i", str(table), "-o", str(tmp_path / "h"), "--mode", "add-header", *extra], tmp_path)
    out = tmp_path / "h" / f"{table.stem}_with_header.csv"
    return r, out


def header_of(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return next(csv.reader(f))


def test_field_dictionary_matches_installed_blast():
    m = load_module()
    help_text = subprocess.run(["blastn", "-help"], env=ENV, capture_output=True, text=True).stdout
    sec = help_text.split("supported format specifiers for options 6, 7 and 10")[1].split("When not provided")[0]
    import re
    codes = re.findall(r"^\s+(\w+) means", sec, flags=re.M)
    assert set(codes) <= set(m.BLAST_FIELDS), set(codes) - set(m.BLAST_FIELDS)


@pytest.mark.parametrize("outfmt,expect_cols", [
    ("6", "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore"),
    ("7 qseqid sseqid pident qcovs stitle evalue", "qseqid sseqid pident qcovs stitle evalue"),
    ("10", "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore"),
    ("20 qseqid sacc nident gaps", "qseqid sacc nident gaps"),
])
def test_add_header_auto_detects(db, tmp_path, outfmt, expect_cols):
    table = blast_table(db, tmp_path, outfmt, "t.out")
    r, out = add_header(table, tmp_path)
    assert r.returncode == 0, r.stdout[-1500:]
    assert header_of(out) == expect_cols.split()
    n_data = sum(1 for x in table.read_text().splitlines() if x and not x.startswith("#"))
    if outfmt.startswith("20"):
        n_data -= 1  # its own header row
    assert len(read_csv_sig(out)) == n_data


def read_csv_sig(p):
    with open(p, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def test_add_header_explicit_outfmt_and_labels(db, tmp_path):
    table = blast_table(db, tmp_path, "6 qseqid sseqid evalue bitscore qcovs", "custom.tsv")
    r, _ = add_header(table, tmp_path)  # 5 headerless columns: cannot be known
    assert r.returncode != 0 and "cannot know the columns" in r.stdout
    r, out = add_header(table, tmp_path, "--outfmt_blast", "6 qseqid sseqid evalue bitscore qcovs",
                        "--header_style", "both")
    assert r.returncode == 0, r.stdout[-1500:]
    assert header_of(out)[2:] == ["evalue (Expect value)", "bitscore (Bit score)",
                                  "qcovs (Query Coverage Per Subject)"]


def test_add_header_rejects_wrong_guess(db, tmp_path):
    # 12 columns but not 'std': guessing std would mislabel -> type check must refuse
    table = blast_table(db, tmp_path, "6 qseqid sseqid stitle sacc saccver qacc qaccver sstrand evalue bitscore qcovs pident", "x.tsv")
    r, out = add_header(table, tmp_path)
    assert r.returncode != 0 and "guess is wrong" in r.stdout
    assert not out.exists()


def test_pipeline_headers_and_topn_resume(db, tmp_path):
    out = tmp_path / "o"
    args = local_args(DATA / "query.fasta", out, db, "--no_smk", "--mode", "blast2topN", "--header_style", "label")
    r = run(args, tmp_path)
    assert r.returncode == 0, r.stdout[-2000:]
    merged = (out / "blast_results" / "query_blast.tsv").read_text().splitlines()
    assert merged[0].split("\t")[:3] == ["qseqid", "sseqid", "stitle"]
    assert header_of(out / "query_results.csv")[:4] == ["Query Seq-id", "Subject Seq-id", "Subject Title",
                                                         "Percentage of identical matches"]
    # resume re-reads the TopN TSV; its header row must not become a data row
    (out / "query_results.csv").unlink()
    (out / "query_completed.txt").unlink()
    r = run(args + ["--resume", "--header_style", "code"], tmp_path)
    assert r.returncode == 0, r.stdout[-2000:]
    rows = read_csv_sig(out / "query_results.csv")
    assert all(row["qseqid"] != "qseqid" for row in rows)
    assert len({row["qseqid"] for row in rows}) == (DATA / "query.fasta").read_text().count(">")
