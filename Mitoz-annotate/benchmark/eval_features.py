"""Per-gene boundary agreement of an annotation with the NCBI record (PCGs and the two rRNAs).

Gene names are normalised (cox1/COX1/COI -> cox1; rrnS/12S/s-rRNA -> rrnS; rrnL/16S/l-rRNA -> rrnL).
For each gene present in both: 5' end same?, 3' end same?, distance. Also counts overlaps of rRNA
with neighbouring features (Fourdrilis criterion vi: rRNA not extended into flanking genes).

Usage: python3 eval_features.py <ncbi.gb> <annotation.gbf or dir> [label]
"""
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mito_names import ALIAS, PCG, PRODUCT  # noqa: E402


def norm(f):
    g = re.sub(r"[^a-z0-9]", "", (f.qualifiers.get("gene", [""])[0] or "").lower())
    g = ALIAS.get(g, g)
    if f.type == "rRNA" and g not in ("rrnS", "rrnL"):
        g = PRODUCT.get(f.qualifiers.get("product", [""])[0].lower(), g)
    return g


def to_ref(rec, ref):
    """Map rec coordinates onto ref (same circular molecule, possibly rotated and/or reverse-
    complemented). Returns f(pos, strand) -> (pos, strand), or None if no mapping found."""
    q, t, L = str(rec.seq).upper(), str(ref.seq).upper(), len(ref.seq)
    if len(q) != L:
        return None
    tt = t + t
    k = tt.find(q[:60])
    if k >= 0:
        return lambda p, st: (((p - 1 + k) % L) + 1, st)
    rc = str(rec.seq.reverse_complement()).upper()
    k = tt.find(rc[:60])
    if k >= 0:
        return lambda p, st: (((L - p + k) % L) + 1, -st)
    return None


def genes(rec, mapf=None):
    out = {}
    for f in rec.features:
        if f.type not in ("CDS", "rRNA"):
            continue
        g = norm(f)
        if g in PCG or g in ("rrnS", "rrnL"):
            s, e = int(f.location.start) + 1, int(f.location.end)
            five, three = (s, e) if f.location.strand == 1 else (e, s)
            st = f.location.strand
            if mapf:
                (five, st), (three, _) = mapf(five, st), mapf(three, f.location.strand)
            out.setdefault(g, (five, three, st, s, e))
    return out


def neighbours_overlap(rec, g, s, e):
    n = 0
    for f in rec.features:
        if f.type in ("tRNA", "CDS", "rRNA") and norm(f) != g:
            fs, fe = int(f.location.start) + 1, int(f.location.end)
            if fs <= e and fe >= s and not (fs <= s and fe >= e):
                n += 1
    return n


ncbi = {r.id.split(".")[0]: r for r in SeqIO.parse(sys.argv[1], "genbank")}
src = Path(sys.argv[2])
files = sorted(src.rglob("*.gbf")) if src.is_dir() else [src]
label = sys.argv[3] if len(sys.argv) > 3 else src.name
stats = defaultdict(Counter)
dist = defaultdict(list)
for f in files:
    for rec in SeqIO.parse(str(f), "genbank"):
        acc = rec.id.split(".")[0]
        if acc not in ncbi:
            continue
        mapf = to_ref(rec, ncbi[acc])
        if mapf is None:
            stats["unmapped"]["n"] += 1
            continue
        a, b = genes(rec, mapf), genes(ncbi[acc])
        for g in sorted(set(a) & set(b)):
            cls = "PCG" if g in PCG else g
            stats[cls]["n"] += 1
            stats[cls]["5' same"] += a[g][0] == b[g][0]
            stats[cls]["3' same"] += a[g][1] == b[g][1]
            dist[cls].append(abs(a[g][0] - b[g][0]) + abs(a[g][1] - b[g][1]))
            if cls != "PCG":
                stats[cls]["overlaps neighbour"] += neighbours_overlap(rec, g, a[g][3], a[g][4]) > 0
print(f"{label}: records not mappable to NCBI sequence: {stats['unmapped']['n']}")
for cls in ("PCG", "rrnS", "rrnL"):
    c = stats[cls]
    if c["n"]:
        print(f"{label:12s} {cls:5s} n={c['n']:4d}  5' same {c[chr(53) + chr(39) + ' same']:4d}  3' same {c[chr(51) + chr(39) + ' same']:4d}"
              f"  median |diff| {median(dist[cls]):.0f} nt" + (f"  overlaps neighbour {c['overlaps neighbour']}" if cls != "PCG" else ""))
