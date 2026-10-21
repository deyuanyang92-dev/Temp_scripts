"""rRNA detection on the dev split: MITOS2 CM vs per-phylum profile HMM (nhmmer) vs per-phylum blastn.

References come from fetch_refs.py (refs/<phylum>.gb, benchmark accessions excluded). For each dev query
the references of its own species are dropped (mode "species") or of its own genus (mode "genus"), so a
query never sees itself. Hits of one gene closer than GAP nt on one strand are merged into one span (this
also merges MITOS2-style fragments); the best-scoring span per gene is the prediction.
Scoring as benchmark/mito_bench_eval.py: a prediction matches an NCBI rRNA of the same name when it covers
>= 50 % of it; ends exact / within 3 / within 10 nt.

Run: python3 -I rr_detect.py <work dir> <species|genus> [phylum,phylum]   (writes <work>/rr_detect_<mode>.tsv and prints a summary)
"""
import csv
import random
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

BIN = "/root/miniconda3/envs/mitos/bin/"
MAFFT = "/root/miniconda3/bin/mafft"
SET = "/mnt/n/yzz-分类工作流/temp_scripts/Mitoz-annotate/benchmark/sets/mito_bench_2026-10.tsv"
TRUTH = Path("/mnt/n/yzz-分类工作流/projects/mito_benchmark/gb")
SCR = Path("/tmp/claude-0/-mnt-n-yzz------/934ca37c-5895-4358-b309-900011b04d11/scratchpad/rrna/cm")
GAP, NSUB, CPU = 300, 60, 6


def rname(f):
    s = " ".join(f.qualifiers.get("product", []) + f.qualifiers.get("gene", [])).lower()
    if any(k in s for k in ("12s", "s-rrna", "small", "rns", "rrns", "srrna", "srna", "ssu")):
        return "rrnS"
    if any(k in s for k in ("16s", "l-rrna", "large", "rnl", "rrnl", "lrrna", "lrna", "lsu")):
        return "rrnL"
    return None


def load_refs(work):
    refs = defaultdict(list)  # (phylum, gene) -> [(organism, genus, seq)]
    for gb in sorted((work / "refs").glob("*.gb")):
        for r in SeqIO.parse(str(gb), "genbank"):
            org = r.annotations.get("organism", "")
            for f in r.features:
                n = rname(f) if f.type == "rRNA" else None
                if not n:
                    continue
                s = str(f.extract(r.seq)).upper()
                if 300 <= len(s) <= 2500 and s.count("N") <= 0.01 * len(s):
                    refs[(gb.stem, n)].append((org, org.split()[0], s))
    return refs


def truth(acc):
    r = SeqIO.read(str(TRUTH / f"{acc}.gb"), "genbank")
    out = []
    for f in r.features:
        if f.type == "rRNA" and rname(f):
            out.append((rname(f), int(f.location.start), int(f.location.end), f.location.strand))
    return out


def cluster(hits):
    """hits: (start, end, strand, score) -> best merged span (start, end, strand) or None."""
    best = None
    for st in (1, -1):
        hs = sorted(h for h in hits if h[2] == st)
        cur = None
        for h in hs:
            if cur and h[0] - cur[1] <= GAP:
                cur = [cur[0], max(cur[1], h[1]), st, cur[3] + h[3]]
            else:
                if cur and (best is None or cur[3] > best[3]):
                    best = cur
                cur = list(h)
        if cur and (best is None or cur[3] > best[3]):
            best = cur
    return tuple(best[:3]) if best else None


def cm_pred(acc, gene, merge):
    tag = "mt12" if gene == "rrnS" else "mt16"
    for mode in ("g", "l"):
        hits = []
        for line in open(SCR / "out" / f"{acc}.{tag}.{mode}.tbl"):
            if line.startswith("#"):
                continue
            c = line.split()
            a, b, st, sc, ev = int(c[7]), int(c[8]), c[9], float(c[14]), float(c[15])
            if ev > 0.01:
                continue
            s, e = (a - 1, b) if st == "+" else (b - 1, a)
            hits.append((s, e, 1 if st == "+" else -1, sc))
        if hits:
            if merge:
                return cluster(hits)
            h = max(hits, key=lambda x: x[3])
            return h[:3]
    return None


def build_db(work, key, seqs, tmp):
    """Write refs FASTA (all) and a profile HMM (from a NSUB subsample); cached by key."""
    d = work / "models" / key
    hmm, fa = d / "m.hmm", d / "r.fa"
    if hmm.exists():
        return hmm, fa
    d.mkdir(parents=True, exist_ok=True)
    fa.write_text("".join(f">r{i}\n{s}\n" for i, s in enumerate(seqs)))
    rnd = random.Random(1)
    sub = seqs if len(seqs) <= NSUB else rnd.sample(seqs, NSUB)
    (d / "sub.fa").write_text("".join(f">r{i}\n{s}\n" for i, s in enumerate(sub)))
    aln = subprocess.run([MAFFT, "--auto", "--thread", str(CPU), "--quiet", str(d / "sub.fa")],
                         capture_output=True, text=True, check=True).stdout
    (d / "a.fa").write_text(aln)
    subprocess.run([BIN + "hmmbuild", "--dna", "--cpu", str(CPU), "-n", key, str(hmm), str(d / "a.fa")],
                   capture_output=True, check=True)
    subprocess.run([BIN + "makeblastdb", "-dbtype", "nucl", "-in", str(fa)], capture_output=True, check=True)
    return hmm, fa


def hmm_pred(hmm, qfa, tmp):
    tbl = tmp / "n.tbl"
    subprocess.run([BIN + "nhmmer", "--dna", "--cpu", str(CPU), "-E", "1e-3", "--tblout", str(tbl), "-o", "/dev/null",
                    str(hmm), str(qfa)], check=True)
    hits = []
    for line in open(tbl):
        if line.startswith("#"):
            continue
        c = line.split()
        a, b, st, sc = int(c[8]), int(c[9]), c[11], float(c[13])  # envelope from/to
        s, e = (a - 1, b) if st == "+" else (b - 1, a)
        hits.append((s, e, 1 if st == "+" else -1, sc))
    return cluster(hits)


def blast_pred(fa, qfa):
    out = subprocess.run([BIN + "blastn", "-task", "blastn", "-evalue", "1e-5", "-num_threads", str(CPU),
                          "-max_target_seqs", "500", "-outfmt", "6 qstart qend sstart send bitscore",
                          "-query", str(qfa), "-db", str(fa)], capture_output=True, text=True, check=True).stdout
    hits = []
    for line in out.splitlines():
        qs, qe, ss, se, bs = line.split()
        st = 1 if int(se) > int(ss) else -1
        hits.append((int(qs) - 1, int(qe), st, float(bs)))
    # one HSP set per reference overlaps heavily; cluster() unions them and scores by summed bitscore
    return cluster(hits)


def score(t, p):
    if p is None:
        return "missed", ""
    n, s, e, st = t
    ov = max(0, min(e, p[1]) - max(s, p[0]))
    if ov < 0.5 * (e - s):
        return "missed", f"{p[0] + 1}..{p[1]}"
    d = max(abs(s - p[0]), abs(e - p[1]))
    return ("exact" if d == 0 else "le3" if d <= 3 else "le10" if d <= 10 else "gt10"), f"{p[0] + 1}..{p[1]}"


def main():
    work, mode = Path(sys.argv[1]), sys.argv[2]
    rows = [r for r in csv.DictReader(open(SET), delimiter="\t") if r["split"] == "dev"
            and (len(sys.argv) < 4 or r["phylum"] in sys.argv[3].split(","))]
    refs = load_refs(work)
    tmp = work / "tmp"
    tmp.mkdir(exist_ok=True)
    out = []
    for r in rows:
        acc, phy = r["accession"], r["phylum"]
        qfa = SCR / "fa" / f"{acc}.fa"
        for gene in ("rrnS", "rrnL"):
            pool = refs.get((phy, gene), [])
            keep = [s for o, g, s in pool if (o != r["species"] if mode == "species" else g != r["genus"])]
            preds = {"cm": cm_pred(acc, gene, False), "cm_merge": cm_pred(acc, gene, True)}
            t0 = time.time()
            if len(keep) >= 3:
                excl = (r["species"] if mode == "species" else r["genus"]).replace(" ", "_")
                hmm, fa = build_db(work, f"{phy}.{gene}.{mode}.{excl}", keep, tmp)
                preds["hmm"] = hmm_pred(hmm, qfa, tmp)
                preds["blast"] = blast_pred(fa, qfa)
            else:
                preds["hmm"] = preds["blast"] = None
            preds["cm+hmm"] = preds["cm_merge"] or preds["hmm"]
            preds["cm+blast"] = preds["cm_merge"] or preds["blast"]
            ts = [t for t in truth(acc) if t[0] == gene]
            for t in ts:
                for m, p in preds.items():
                    c, pos = score(t, p)
                    out.append(dict(acc=acc, phylum=phy, gene=gene, nref=len(keep), method=m, truth=f"{t[1] + 1}..{t[2]}",
                                    pred=pos, cls=c))
            print(acc, phy, gene, len(keep), f"{time.time() - t0:.0f}s", {m: score(ts[0], p)[0] for m, p in preds.items()} if ts else "-",
                  flush=True)
    f = work / f"rr_detect_{mode}.tsv"
    with open(f, "w") as h:
        w = csv.DictWriter(h, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    tot = defaultdict(lambda: defaultdict(int))
    for o in out:
        tot[o["method"]][o["cls"]] += 1
    for m, d in tot.items():
        n = sum(d.values())
        print(f"{m:9s} found {n - d['missed']}/{n}  exact {d['exact']}  <=3 {d['exact'] + d['le3']}  <=10 {d['exact'] + d['le3'] + d['le10']}")


if __name__ == "__main__":
    main()
