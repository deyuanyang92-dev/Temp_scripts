"""Real-database benchmark: cost of one database pass vs. per-query cost, thread mode, process split.

python benchmarks/benchmark_nt_passes.py --db-dir /path/to/blast_db/nt --volumes 10 --threads 8

Uses the first N volumes of nt (warm them into the page cache once) and mutated copies of
tests/data/query.fasta (3% substitutions, so every query is distinct). All configurations must
return the same HSP set; only wall time differs. Results in docs/benchmark-nt-passes.jsonl.
"""
import argparse
import concurrent.futures as cf
import hashlib
import json
import random
import subprocess
import tempfile
import time
from pathlib import Path

FMT = "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
DATA = Path(__file__).resolve().parents[1] / "tests" / "data" / "query.fasta"


def make_queries(root, n):
    random.seed(1)
    seqs = [s.split("\n", 1)[1].replace("\n", "") for s in DATA.read_text().split(">")[1:]]
    path = root / f"q{n}.fa"
    with open(path, "w") as f:
        for i in range(n):
            s = "".join(random.choice("ACGT") if random.random() < 0.03 else c for c in seqs[i % len(seqs)])
            f.write(f">q{i}\n{s}\n")
    return path


def split(src, n):
    recs = src.read_text().split(">")[1:]
    k = -(-len(recs) // n)
    out = []
    for i in range(0, len(recs), k):
        p = src.with_name(f"{src.stem}_{n}_{i}.fa")
        p.write_text("".join(">" + r for r in recs[i:i + k]))
        out.append(p)
    return out


def blast(q, db, threads, mt):
    out = q.with_suffix(f".t{threads}m{mt}.tsv")
    subprocess.run(["blastn", "-query", str(q), "-db", db, "-outfmt", FMT, "-out", str(out),
                    "-evalue", "1e-5", "-max_target_seqs", "100", "-num_threads", str(threads),
                    "-mt_mode", str(mt)], check=True, capture_output=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db-dir", required=True)
    ap.add_argument("--volumes", type=int, default=10)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    db = " ".join(str(Path(a.db_dir) / f"nt.{i:03d}") for i in range(a.volumes))
    root = Path(tempfile.mkdtemp(prefix="b2m-nt-bench-"))
    q = {n: make_queries(root, n) for n in (6, 96, 960)}
    blast(q[6], db, a.threads, 0)  # warm page cache
    t = a.threads
    # (name, queries, batches, processes, mt_mode)
    for name, n, nb, procs, mt in [
        ("6 queries, 1 proc", 6, 1, 1, 0), ("96 queries, 1 proc", 96, 1, 1, 0),
        ("960 queries, 1 batch, 1 proc", 960, 1, 1, 0), ("960 queries, 1 batch, mt_mode 1", 960, 1, 1, 1),
        ("960 queries, 2 batches, 2 procs", 960, 2, 2, 0), ("960 queries, 20 batches, 4 procs", 960, 20, 4, 0),
    ]:
        qs = split(q[n], nb)
        start = time.perf_counter()
        with cf.ThreadPoolExecutor(procs) as ex:
            outs = list(ex.map(lambda x: blast(x, db, t // procs, mt), qs))
        rows = sorted(line for o in outs for line in o.read_text().splitlines())
        print(json.dumps(dict(name=name, threads=t, seconds=round(time.perf_counter() - start, 1), rows=len(rows),
                              sha1=hashlib.sha1("\n".join(rows).encode()).hexdigest()[:12])), flush=True)


if __name__ == "__main__":
    main()
