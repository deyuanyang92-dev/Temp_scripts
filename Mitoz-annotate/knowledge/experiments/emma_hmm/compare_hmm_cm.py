"""Compare result_emma_hmm_dev.txt (EMMA vertebrate nhmmer HMMs) with result_mitos2_cm_dev.txt (MITOS2 RefSeq-89 CMs, cmsearch -g)
on the same 93 dev genomes.  Keys (accession, gene, truth length) collapse 2 duplicated features, so 184 of 186 are compared.
Usage: python3 -I compare_hmm_cm.py > compare_hmm_cm.txt
"""
import collections
TARGET = {'Chaetognatha', 'Ctenophora', 'Placozoa', 'Xenacoelomorpha', 'Kinorhyncha', 'Gnathostomulida', 'Nematoda', 'Tardigrada'}


def load(fn):
    d, on = {}, False
    for line in open(fn):
        if line.startswith("## per-record"):
            on = True
            continue
        if on:
            b = line.rstrip("\n").split("\t")
            d[(b[1], b[2], b[3])] = (b[0], float(b[5]), float(b[6]), float(b[7]))
    return d


h = load("result_emma_hmm_dev.txt")
c = load("result_mitos2_cm_dev.txt")
print("keys", len(h), len(c), "(features collapsed by duplicate keys:", 186 - len(h), ")")
for name, i in (("E<1e-10", 1), ("E<1e-5", 2), ("E<10", 3)):
    for label, sel in (("all", lambda p: True), ("target8", lambda p: p in TARGET)):
        ks = [k for k in h if sel(h[k][0])]
        nh = sum(h[k][i] >= 0.5 for k in ks)
        nc = sum(c[k][i] >= 0.5 for k in ks)
        both = sum(h[k][i] >= 0.5 and c[k][i] >= 0.5 for k in ks)
        honly = sum(h[k][i] >= 0.5 and not c[k][i] >= 0.5 for k in ks)
        conly = sum(c[k][i] >= 0.5 and not h[k][i] >= 0.5 for k in ks)
        print(f"{name} {label}: n={len(ks)} HMM={nh} CM={nc} both={both} HMM_only={honly} CM_only={conly}")
print("target-8 features missed by the CMs at E<1e-10, by phylum:")
miss = collections.Counter(h[k][0] for k in h if h[k][0] in TARGET and c[k][1] < 0.5)
print(dict(miss))
