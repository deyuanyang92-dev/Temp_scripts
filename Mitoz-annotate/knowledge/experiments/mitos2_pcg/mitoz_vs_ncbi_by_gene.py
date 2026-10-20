"""Per-gene boundary agreement of MitoZ GenBank output with NCBI records (knowledge experiment, 2026-10-20).

Usage: python3 mitoz_vs_ncbi_by_gene.py ncbi.gb 'mitoz_dir/*/*.gbf' [junctions.tsv] > by_gene.tsv
Optional third argument: write the nad4l/nad4 overlap (nt, + strand, adjacent genes only) in the NCBI record and in the MitoZ output.
MitoZ was run on the NCBI sequences, but the sequence may be rotated or reverse-complemented; the MitoZ sequence is located in
the NCBI sequence (first 60 nt, forward or reverse complement) and its coordinates are converted. Single-interval NCBI CDS only.
Columns: gene, ncbi_cds, mitoz_found, wrong_locus, five_prime_identical, three_prime_identical, both_identical.
"""
import collections
import glob
import sys

from Bio import SeqIO
from Bio.Seq import Seq

from pcg_names import pcg_name as gname

GENES = ["atp8", "atp6", "nad4l", "nad4", "nad6", "cox1", "cox2", "cox3", "cob", "nad1", "nad2", "nad3", "nad5"]


ncbi = {r.id.split(".")[0]: r for r in SeqIO.parse(sys.argv[1], "genbank")}
mz = {}
for f in glob.glob(sys.argv[2]):
    for r in SeqIO.parse(f, "genbank"):
        mz[r.id.split(".")[0]] = r
res = collections.defaultdict(collections.Counter)
used = 0
junc = []
for acc, rm in mz.items():
    rn = ncbi[acc]
    cn = {gname(f): (int(f.location.start) + 1, int(f.location.end), f.location.strand)
          for f in rn.features if f.type == "CDS" and len(f.location.parts) == 1 and gname(f)}
    if not cn:
        continue
    sm, sn = str(rm.seq).upper(), str(rn.seq).upper()
    L = len(sn)
    if len(sm) != L:
        continue
    k = (sn + sn).find(sm[:60])
    flip = False
    if k < 0:
        k = (sn + sn).find(str(Seq(sm).reverse_complement())[:60])
        flip = True
        if k < 0:
            continue
    used += 1
    cm = collections.defaultdict(list)
    for f in rm.features:
        if f.type == "CDS" and gname(f):
            s0, e0 = int(str(f.location.start).lstrip("<>")) + 1, int(str(f.location.end).lstrip("<>"))
            st = f.location.strand
            if flip:
                s0, e0, st = L - e0 + 1, L - s0 + 1, -st
            s1, e1 = (s0 - 1 + k) % L + 1, (e0 - 1 + k) % L + 1
            if s1 <= e1:
                cm[gname(f)].append((s1, e1, st))
    if "nad4l" in cn and "nad4" in cn and cm.get("nad4l") and cm.get("nad4") and cn["nad4l"][2] == cn["nad4"][2] == 1:
        ov_n = cn["nad4l"][1] - cn["nad4"][0] + 1
        ov_m = cm["nad4l"][0][1] - cm["nad4"][0][0] + 1
        if abs(ov_n) < 100 and abs(ov_m) < 100:
            junc.append((acc, ov_n, ov_m))
    for g, a in cn.items():
        c = res[g]
        c["ncbi_cds"] += 1
        b = cm.get(g)
        if not b:
            continue
        b = max(b, key=lambda x: min(x[1], a[1]) - max(x[0], a[0]))
        if min(b[1], a[1]) - max(b[0], a[0]) < 0 or b[2] != a[2]:
            c["wrong_locus"] += 1
            continue
        c["mitoz_found"] += 1
        d5 = b[0] - a[0] if a[2] == 1 else a[1] - b[1]
        d3 = b[1] - a[1] if a[2] == 1 else a[0] - b[0]
        c["five_prime_identical"] += d5 == 0
        c["three_prime_identical"] += d3 == 0
        c["both_identical"] += d5 == 0 and d3 == 0
print(f"# records compared: {used}")
print("gene\tncbi_cds\tmitoz_found\twrong_locus\tfive_prime_identical\tthree_prime_identical\tboth_identical")
for g in GENES:
    c = res[g]
    print("\t".join([g] + [str(c[k]) for k in ("ncbi_cds", "mitoz_found", "wrong_locus", "five_prime_identical",
                                               "three_prime_identical", "both_identical")]))

if len(sys.argv) > 3:
    with open(sys.argv[3], "w") as out:
        out.write("acc\tncbi_overlap_nt\tmitoz_overlap_nt\n")
        for row in junc:
            out.write("\t".join(map(str, row)) + "\n")
