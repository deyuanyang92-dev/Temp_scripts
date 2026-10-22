"""Online (NCBI qblast) benchmark: old per-query submission vs. new one-request-per-batch.

python benchmarks/benchmark_online_qblast.py --email you@example.org [--n 10] [--old-n 3]

Sends real requests to NCBI; keep --n small and run outside US peak hours. "old" replays the
pre-2026-10-23 loop (one qblast per sequence + 12 s sleep, 3.1 s limiter) on the first --old-n
queries; "new" runs the package's run_biopython_qblast_one_batch on all --n queries.
Also probes MEGABLAST=on vs. the URL-API default for the same batch.
"""
import argparse
import json
import random
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from Bio.Blast import NCBIWWW, NCBIXML  # noqa: E402

from blast2metadata import cli, engines, ratelimit  # noqa: E402
from blast2metadata.fastautil import BatchInfo  # noqa: E402


def make_queries(path, n):
    random.seed(1)
    seqs = [s.split("\n", 1)[1].replace("\n", "") for s in (ROOT / "tests/data/query.fasta").read_text().split(">")[1:]]
    with open(path, "w") as f:
        for i in range(n):
            s = "".join(random.choice("ACGT") if random.random() < 0.03 else c for c in seqs[i % len(seqs)])
            f.write(f">q{i}\n{s}\n")


def old_loop(fasta, n, args):
    """Faithful replay of the previous per-sequence code path (timing only)."""
    from Bio import SeqIO
    pairs = set()
    for rec in list(SeqIO.parse(str(fasta), "fasta"))[:n]:
        time.sleep(0)  # old limiter (3.1 s) is always satisfied because of the 12 s sleep below
        h = NCBIWWW.qblast(program="blastn", database="nt", sequence=str(rec.seq),
                           expect=args.evalue, hitlist_size=args.max_target_seqs)
        for br in NCBIXML.parse(h):
            pairs |= {(rec.id, a.accession) for a in br.alignments}
        h.close()
        time.sleep(12)
    return pairs


def new_batch(fasta, tmp, args, label):
    out = tmp / label
    _, ok, path, cols = engines.run_biopython_qblast_one_batch(BatchInfo(label, fasta), "nucleotide", args, out,
                                                               engines.BLAST_HEADERS)
    assert ok
    rows = [line.split("\t") for line in path.read_text().splitlines()]
    return {(r[0], r[1]) for r in rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True)
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--old-n", type=int, default=3)
    ap.add_argument("--max-target-seqs", type=int, default=100)
    a = ap.parse_args()
    tmp = Path(tempfile.mkdtemp(prefix="b2m-online-"))
    fasta = tmp / "q.fa"
    make_queries(fasta, a.n)
    base = ["-i", str(fasta), "--no_smk", "-bm", "online", "--online_engine", "biopython",
            "--email", a.email, "--max_target_seqs", str(a.max_target_seqs)]
    args = cli.parse_arguments(base)
    engines.entrez_setup(args)
    ratelimit.BLAST_LIMITER = ratelimit._make_blast_limiter()

    def report(name, n, seconds, pairs):
        print(json.dumps(dict(name=name, queries=n, seconds=round(seconds, 1),
                              per_query=round(seconds / n, 1), pairs=len(pairs))), flush=True)
        return pairs

    t = time.time(); new = new_batch(fasta, tmp, args, "new_default")
    report("new: 1 request per batch (URL-API default task)", a.n, time.time() - t, new)
    t = time.time(); mb = new_batch(fasta, tmp, cli.parse_arguments(base + ["--extra_para2_blast", "-megablast on"]),
                                    "new_megablast")
    report("new: 1 request per batch, MEGABLAST=on", a.n, time.time() - t, mb)
    t = time.time(); old = old_loop(fasta, a.old_n, args)
    report("old: 1 request per query + 12 s sleep", a.old_n, time.time() - t, old)
    ids = {f"q{i}" for i in range(a.old_n)}
    same = {p for p in new if p[0] in ids} == old
    print(json.dumps(dict(old_vs_new_same_pairs_on_shared_queries=same,
                          default_vs_megablast_same_pairs=new == mb)), flush=True)


if __name__ == "__main__":
    main()
