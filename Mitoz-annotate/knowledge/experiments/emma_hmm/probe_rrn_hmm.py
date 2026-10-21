"""Probe: how well do EMMA's vertebrate-built nhmmer rRNA HMMs (emma-models rrn/all_rrn.hmm) find the NCBI rRNAs
of the 93 DEV records of benchmark/sets/mito_bench_2026-10.tsv?  (test split not used.)
Mimics only the search step of EMMA (genome + first 100 nt, plain `nhmmer --tblout`, E < 1e-10 filter), NOT its
merge/snapping code.
Usage HMM mode: python3 -I probe_rrn_hmm.py hmm <all_rrn.hmm> <workdir> <nhmmer_bin> > result.txt
Usage CM mode (comparison with the MITOS2 CMs, glocal like MITOS2's first pass):
       python3 -I probe_rrn_hmm.py cm <rrnS.cm> <rrnL.cm> <workdir> <cmsearch_bin> > result.txt
"""
import csv, collections, subprocess, sys, os, statistics
from Bio import SeqIO

MODE = sys.argv[1]
if MODE == "hmm":
    HMM, WORK, NHMMER = sys.argv[2:5]
else:
    CMS, WORK, CMSEARCH = {"rrnS": sys.argv[2], "rrnL": sys.argv[3]}, sys.argv[4], sys.argv[5]
R = "/mnt/n/yzz-分类工作流/projects/mito_benchmark"
BENCH = "/mnt/n/yzz-分类工作流/temp_scripts/Mitoz-annotate/benchmark/sets/mito_bench_2026-10.tsv"
rows = [r for r in csv.DictReader(open(BENCH), delimiter="\t") if r["split"] == "dev"]
os.makedirs(WORK, exist_ok=True)


def which(f):
    nm = (f.qualifiers.get("gene") or f.qualifiers.get("product") or ["?"])[0]
    prod = (f.qualifiers.get("product") or [""])[0].lower()
    gn = nm.lower()
    if "12s" in prod or "small" in prod or "rrns" in gn or "12s" in gn or "rns" in gn or "s-rrna" in gn:
        return "rrnS"
    if "16s" in prod or "large" in prod or "rrnl" in gn or "16s" in gn or "rnl" in gn or "l-rrna" in gn:
        return "rrnL"
    return "?"


MODEL = {"12srna": "rrnS", "16srna": "rrnL"}
res = []   # (phylum, acc, gene, truth_len, best_E, cov_union_1e10, cov_union_1e5, cov_union_E10, nhits_1e10)
for r in rows:
    acc = r["accession"]
    try:
        rec = SeqIO.read(f"{R}/gb/{acc}.gb", "genbank")
    except Exception as e:
        print("skip", acc, e, file=sys.stderr)
        continue
    L = len(rec.seq)
    seq = str(rec.seq).upper()
    fa = f"{WORK}/{acc}.ext.fa"
    open(fa, "w").write(f">{acc}\n{seq + seq[:100]}\n")
    hits = []
    if MODE == "hmm":
        tbl = f"{WORK}/{acc}.tbl"
        subprocess.run([NHMMER, "--cpu", "1", "--tblout", tbl, "-o", "/dev/null", HMM, fa], check=True)
        for line in open(tbl):
            if line.startswith("#"):
                continue
            b = line.split()
            a, z = int(b[6]), int(b[7])
            hits.append((MODEL[b[2]], b[11], min(a, z), max(a, z), float(b[12]), float(b[13]), int(b[4]), int(b[5])))
    else:
        for gname, cm in CMS.items():
            tbl = f"{WORK}/{acc}.{gname}.tbl"
            subprocess.run([CMSEARCH, "-g", "--cpu", "1", "--tblout", tbl, "-o", "/dev/null", cm, fa], check=True)
            for line in open(tbl):
                if line.startswith("#"):
                    continue
                b = line.split()
                a, z = int(b[7]), int(b[8])
                hits.append((gname, b[9], min(a, z), max(a, z), float(b[15]), float(b[14]), int(b[5]), int(b[6])))
    for f in rec.features:
        if f.type != "rRNA":
            continue
        g = which(f)
        if g == "?":
            continue
        pos = set()
        for p in f.location.parts:
            pos.update(range(int(p.start) + 1, int(p.end) + 1))
        strand = "+" if f.location.strand != -1 else "-"
        covs = {}
        for thr in (1e-10, 1e-5, 10.0):
            cov = set()
            for (m, st, s, e, E, sc, hf, ht) in hits:
                if m == g and st == strand and E < thr:
                    cov.update(((i - 1) % L) + 1 for i in range(s, e + 1))
            covs[thr] = len(cov & pos) / len(pos)
        bestE = min([h[4] for h in hits if h[0] == g and h[1] == strand] or [float("nan")])
        wrong = [h for h in hits if h[0] != g and h[4] < 1e-10]
        res.append((r["phylum"], acc, g, len(pos), bestE, covs[1e-10], covs[1e-5], covs[10.0],
                    sum(1 for h in hits if h[0] == g and h[1] == strand and h[4] < 1e-10), len(wrong)))

print("# EMMA rRNA HMM probe, dev split; coverage = fraction of the NCBI rRNA (same strand) covered by union of nhmmer hits of the matching model")
print("# cols: phylum n_truth detected(cov>=0.5 at E<1e-10) detected(E<1e-5) detected(E<10) median_cov(E<1e-10) mean_nhits")
bp = collections.defaultdict(list)
for x in res:
    bp[x[0]].append(x)
order = sorted(bp, key=lambda p: -len(bp[p]))
for gene in ("rrnS", "rrnL", "both"):
    print(f"## {gene}")
    for p in order:
        xs = [x for x in bp[p] if gene == "both" or x[2] == gene]
        if not xs:
            continue
        d10 = sum(x[5] >= 0.5 for x in xs)
        d5 = sum(x[6] >= 0.5 for x in xs)
        d0 = sum(x[7] >= 0.5 for x in xs)
        print(f"{p}\t{len(xs)}\t{d10}\t{d5}\t{d0}\t{statistics.median(x[5] for x in xs):.2f}\t{statistics.mean(x[8] for x in xs):.1f}")
    xs = [x for x in res if gene == "both" or x[2] == gene]
    print(f"ALL\t{len(xs)}\t{sum(x[5]>=0.5 for x in xs)}\t{sum(x[6]>=0.5 for x in xs)}\t{sum(x[7]>=0.5 for x in xs)}\t{statistics.median(x[5] for x in xs):.2f}")
print("## per-record (phylum acc gene truth_len bestE cov1e-10 cov1e-5 covE10 nhits1e-10 n_wrong_model_hits)")
for x in sorted(res):
    print("\t".join(f"{v:.3g}" if isinstance(v, float) else str(v) for v in x))
