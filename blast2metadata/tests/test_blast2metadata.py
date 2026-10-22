"""Regression tests for blast2metadata with a real local BLAST+ database.

Run: python -m pytest -q tests/
Needs blastn and makeblastdb on PATH (or BLAST_BIN=<dir>).
"""
import csv
import io
import os
import shutil
import subprocess
import sys
import types
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
SRC = HERE.parent / "src"
SCRIPT = HERE.parent / "run_blast2metadata.py"
sys.path.insert(0, str(SRC))
from blast2metadata import accession, blasttab, cli, engines, fields, ratelimit  # noqa: E402
BLAST_BIN = os.environ.get("BLAST_BIN") or str(Path(sys.executable).parent)
ENV = {**os.environ, "PATH": BLAST_BIN + os.pathsep + os.environ.get("PATH", "")}

pytestmark = pytest.mark.skipif(
    shutil.which("blastn", path=ENV["PATH"]) is None, reason="BLAST+ not installed"
)


def load_module():
    """Flat view of the package API (the tests predate the package split)."""
    m = types.SimpleNamespace()
    for mod in (fields, accession, blasttab, engines, ratelimit, cli):
        m.__dict__.update({k: v for k, v in vars(mod).items() if not k.startswith("__")})
    return m


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
    ratelimit.BLAST_LIMITER = m._RateLimiter(1000, 1.0)
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
    if outfmt.startswith("20") and "20" not in blastn_outfmts():
        pytest.skip("outfmt 20 needs BLAST+ >= 2.13")
    table = blast_table(db, tmp_path, outfmt, "t.out")
    r, out = add_header(table, tmp_path)
    assert r.returncode == 0, r.stdout[-1500:]
    assert header_of(out) == expect_cols.split()
    n_data = sum(1 for x in table.read_text().splitlines() if x and not x.startswith("#"))
    if outfmt.startswith("20"):
        n_data -= 1  # its own header row
    assert len(read_csv_sig(out)) == n_data


def blastn_outfmts():
    h = subprocess.run(["blastn", "-help"], env=ENV, capture_output=True, text=True).stdout
    import re
    return set(re.findall(r"^\s+(\d+) = ", h, flags=re.M))


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


# ── TopN ranking options ─────────────────────────────────────────────────

HITS = [  # one query; values chosen so each criterion picks a different winner
    {"qseqid": "q", "sseqid": "short_perfect", "pident": 100.0, "qcovs": 20, "evalue": 1e-30, "bitscore": 200.0},
    {"qseqid": "q", "sseqid": "full_99",       "pident": 99.0,  "qcovs": 100, "evalue": 0.0,  "bitscore": 1090.0},
    {"qseqid": "q", "sseqid": "full_cov_low",  "pident": 90.0,  "qcovs": 100, "evalue": 0.0,  "bitscore": 1100.0},
    {"qseqid": "q", "sseqid": "mid",           "pident": 98.0,  "qcovs": 95,  "evalue": 0.0,  "bitscore": 1000.0},
]


def ranked(m, spec, **kw):
    return [r["sseqid"] for r in m.get_topN_hits_py(HITS, 10, spec, **kw)]


def test_sort_presets_and_custom_columns():
    m = load_module()
    assert ranked(m, "bitscore") == ["full_cov_low", "full_99", "mid", "short_perfect"]
    assert ranked(m, "identity") == ["short_perfect", "full_99", "mid", "full_cov_low"]
    assert ranked(m, "coverage") == ["full_99", "full_cov_low", "mid", "short_perfect"]
    # combined: % of query identical = pident*qcovs/100 -> 99, 93.1, 90, 20
    assert ranked(m, "pident_x_qcovs") == ["full_99", "mid", "full_cov_low", "short_perfect"]
    top = m.get_topN_hits_py(HITS, 1, "pident_x_qcovs")[0]
    assert top["pident_x_qcovs"] == 99.0
    # custom order; aliases; evalue ascending
    assert ranked(m, "qcovs,pident") == ranked(m, "coverage,identity") == ["full_99", "full_cov_low", "mid", "short_perfect"]
    assert ranked(m, "evalue,bitscore")[:3] == ["full_cov_low", "full_99", "mid"]


def test_sort_filters():
    m = load_module()
    assert ranked(m, "identity", min_qcovs=90) == ["full_99", "mid", "full_cov_low"]
    assert ranked(m, "bitscore", min_pident=97, min_qcovs=96) == ["full_99"]


def test_sort_spec_rejects_unknown_column(tmp_path):
    r = run(["-i", "x.fasta", "--no_smk", "--sort_by", "pident,foo"], tmp_path)
    assert r.returncode == 2 and "unknown --sort_by column" in r.stderr


def test_custom_sort_end_to_end(db, tmp_path):
    out = tmp_path / "o"
    r = run(local_args(DATA / "query.fasta", out, db, "--no_smk", "--mode", "blast2topN",
                       "--sort_by", "pident_x_qcovs", "--min_qcovs", "50", "--topN", "2"), tmp_path)
    assert r.returncode == 0, r.stdout[-2000:]
    rows = read_csv(out / "query_results.csv")
    assert rows and all(float(x["qcovs"]) >= 50 for x in rows)
    scores = {}
    for x in rows:
        scores.setdefault(x["qseqid"], []).append(float(x["pident_x_qcovs"]))
    assert all(v == sorted(v, reverse=True) for v in scores.values())


# ── v3.0: package, engines, gzip, Entrez robustness ─────────────────────

DIAMOND_EXE = os.environ.get("DIAMOND_EXE") or shutil.which("diamond")
MMSEQS_EXE = os.environ.get("MMSEQS_EXE") or shutil.which("mmseqs")


def results_by_query(path):
    out = {}
    for row in read_csv(path):
        out.setdefault(row["qseqid"], []).append(row)
    return out


@pytest.fixture(scope="module")
def protein_db_fasta(tmp_path_factory):
    """COX1 proteins translated from the query genes (invertebrate mito code) -> DIAMOND blastx target."""
    from Bio import SeqIO
    d = tmp_path_factory.mktemp("prot")
    out = d / "cox1.faa"
    with open(out, "w") as fh:
        for r in SeqIO.parse(str(DATA / "query.fasta"), "fasta"):
            best = min((str(r.seq[f:len(r.seq) - (len(r.seq) - f) % 3].translate(table=5)) for f in range(3)),
                       key=lambda p: p.count("*"))
            fh.write(f">{r.id}_prot cytochrome c oxidase subunit I [test]\n{best.rstrip('*')}\n")
    return out


@pytest.mark.skipif(DIAMOND_EXE is None, reason="DIAMOND not installed (set DIAMOND_EXE)")
def test_diamond_engine_blastx_from_fasta_db(protein_db_fasta, tmp_path):
    out = tmp_path / "o"
    # the proteins were translated with the invertebrate mito code -> tell DIAMOND via extra params
    base = ["-i", str(DATA / "query.fasta"), "-o", str(out), "-bm", "local", "--engine", "diamond",
            "--blast_db", str(protein_db_fasta), "--diamond_exe", DIAMOND_EXE, "--no_smk",
            "--mode", "blast2topN", "--topN", "1", "--no_ncbi_metadata", "--query_gencode", "5"]
    r = run(base, tmp_path)
    assert r.returncode == 0, r.stdout[-2000:]
    assert (out / "_engine_db" / "cox1.dmnd").exists()  # built once from the FASTA
    assert header_of(out / "query_results.csv") == fields.BLAST_HEADERS
    by_q = results_by_query(out / "query_results.csv")
    for q, rows in by_q.items():
        assert rows[0]["sseqid"] == q + "_prot", (q, rows)  # own translation is the best hit
        # the frameshifted "_indel" query cannot be translated full length
        assert float(rows[0]["qcovs"]) >= (80 if q.endswith("_indel") else 95), (q, rows)
        assert float(rows[0]["pident"]) >= 99, (q, rows)
    assert "NC_006902.1_COX1" in by_q
    # resume reuses the DB and the batch outputs
    r = run(base + ["--resume"], tmp_path)
    assert r.returncode == 0 and "building database" not in r.stdout


@pytest.mark.skipif(MMSEQS_EXE is None, reason="MMseqs2 not installed (set MMSEQS_EXE)")
def test_mmseqs_engine_nucleotide_matches_blast(db, tmp_path):
    out = tmp_path / "o"
    r = run(["-i", str(DATA / "query.fasta"), "-o", str(out), "-bm", "local", "--engine", "mmseqs",
             "--blast_db", str(DATA / "mito_db.fasta"), "--mmseqs_exe", MMSEQS_EXE, "--no_smk",
             "--mode", "blast2topN", "--topN", "1", "--no_ncbi_metadata"], tmp_path)
    assert r.returncode == 0, r.stdout[-2000:]
    m = load_module()
    assert m.mmseqs_db_type(str(out / "_engine_db" / "mito_db_mmseqsdb")) == "nucleotide"
    assert header_of(out / "query_results.csv") == fields.BLAST_HEADERS
    by_q = results_by_query(out / "query_results.csv")
    assert set(by_q) == {r.split()[0][1:] for r in (DATA / "query.fasta").read_text().splitlines() if r.startswith(">")}
    for q, rows in by_q.items():
        assert q.split("_COX")[0] == rows[0]["sseqid"], (q, rows)  # source genome is the best hit
        assert float(rows[0]["pident"]) >= 98, (q, rows)
        # MMseqs2 reports one alignment per query-target pair, so LC670763 COX1 (spans the genome
        # origin: 2 HSPs in BLAST) only gets ~61% coverage here; the others are full length.
        assert float(rows[0]["qcovs"]) >= (60 if q.startswith("LC670763") else 99), (q, rows)
    # a reverse-strand query must still be found (--strand 2)
    rc = tmp_path / "rc.fasta"
    from Bio import SeqIO
    rec = next(SeqIO.parse(str(DATA / "query.fasta"), "fasta"))
    rc.write_text(f">{rec.id}_rc\n{rec.seq.reverse_complement()}\n")
    r = run(["-i", str(rc), "-o", str(tmp_path / "o2"), "-bm", "local", "--engine", "mmseqs",
             "--blast_db", str(out / "_engine_db" / "mito_db_mmseqsdb"), "--mmseqs_exe", MMSEQS_EXE,
             "--no_smk", "--mode", "blast2topN", "--topN", "1", "--no_ncbi_metadata"], tmp_path)
    assert r.returncode == 0, r.stdout[-2000:]
    rows = read_csv(tmp_path / "o2" / "rc_results.csv")
    assert rows and rows[0]["sseqid"] == rec.id.split("_COX")[0]


def test_engine_requires_local_mode(tmp_path):
    r = run(["-i", str(DATA / "query.fasta"), "-o", str(tmp_path / "o"), "--engine", "mmseqs", "--no_smk"], tmp_path)
    assert r.returncode != 0 and "only works with -bm local" in r.stdout


def test_qcovs_from_spans_matches_blast_rounding():
    m = load_module()
    assert m._qcovs_from_spans([(1, 50), (40, 100)], 100) == 100
    assert m._qcovs_from_spans([(1, 50), (52, 100)], 100) == 99   # 99 covered -> 99
    assert m._qcovs_from_spans([(1, 50), (51, 100)], 100) == 100  # adjacent spans merge
    assert m._qcovs_from_spans([(1, 999)], 1000) == 99            # 99.9% is not 100
    assert m._qcovs_from_spans([(1, 5)], 1000) == 1               # 0.5% rounds half up


def test_gzipped_fasta_input_and_output_stem(db, tmp_path):
    import gzip
    d = tmp_path / "in"
    d.mkdir()
    with gzip.open(d / "sample1.fasta.gz", "wt") as fh:
        fh.write((DATA / "query.fasta").read_text())
    r = run(local_args(d, tmp_path / "o", db, "--no_smk", "--mode", "blast2topN", "--topN", "1"), tmp_path)
    assert r.returncode == 0, r.stdout[-2000:]
    rows = read_csv(tmp_path / "o" / "sample1_results.csv")  # stem without .fasta.gz
    assert len(rows) == (DATA / "query.fasta").read_text().count(">")


def test_duplicate_query_ids_are_warned(db, tmp_path):
    fa = tmp_path / "dup.fasta"
    fa.write_text((DATA / "query.fasta").read_text() * 2)
    r = run(local_args(fa, tmp_path / "o", db, "--no_smk", "--mode", "blast-only"), tmp_path)
    assert r.returncode == 0 and "duplicate sequence ID" in r.stdout


def test_entrez_batch_failure_is_split_so_good_accessions_survive():
    from blast2metadata import metadata
    m = load_module()
    calls = []

    class Handle(list):
        def close(self):
            pass

    class FakeEntrez:
        email = api_key = tool = None

        @staticmethod
        def efetch(db, id, rettype, retmode):
            ids = id.split(",")
            calls.append((db, ids))
            if "BAD1" in ids:
                raise RuntimeError("HTTP Error 400: Bad Request")
            return Handle(ids)

        @staticmethod
        def read(handle):
            return [{"GBSeq_primary-accession": a, "GBSeq_accession-version": a + ".1",
                     "GBSeq_organism": "Org " + a} for a in handle]

    real = metadata.Entrez
    metadata.Entrez = FakeEntrez
    ratelimit.ENTREZ_LIMITER = m._RateLimiter(1000, 1.0)
    try:
        args = m.parse_arguments(["-i", "x", "--no_smk", "--retries", "1", "--entrez_delay_sec", "0"])
        got = metadata.fetch_metadata_batch(["OK1", "OK2", "BAD1", "OK3"], args, batch_size=20)
    finally:
        metadata.Entrez = real
    assert set(got) == {"OK1", "OK2", "OK3"}
    assert got["OK3"]["organism"] == "Org OK3"
    assert ("nucleotide", ["OK1", "OK2", "BAD1", "OK3"]) == calls[0]      # whole batch first
    assert ("nucleotide", ["OK3"]) in calls                                # then split down to singles
    assert any(db == "protein" and ids == ["BAD1"] for db, ids in calls)  # alt DB tried for the bad one


def test_console_entry_points(tmp_path):
    r = run(["--version"], tmp_path)
    assert r.returncode == 0 and r.stdout.strip().startswith("blast2metadata 3.")
    r = run(["--doctor"], tmp_path)
    assert r.returncode == 0 and "blastn" in r.stdout and "biopython" in r.stdout
    if DIAMOND_EXE:
        r = run(["--doctor", "--diamond_exe", DIAMOND_EXE], tmp_path)
        assert DIAMOND_EXE in r.stdout and "diamond      NOT FOUND" not in r.stdout
    r = subprocess.run([sys.executable, "-m", "blast2metadata", "--version"], cwd=tmp_path,
                       env={**ENV, "PYTHONPATH": str(SRC)}, capture_output=True, text=True)
    assert r.returncode == 0 and "blast2metadata 3." in r.stdout
    r = run(["--topN", "0", "-i", "x"], tmp_path)
    assert r.returncode == 2 and "must be >= 1" in r.stderr
