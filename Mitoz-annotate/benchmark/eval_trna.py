"""Evaluate mito_trna.py against NCBI tRNA annotation: each tool, and the selected set at several
min-support levels (independent method families).

Matching: same strand, overlap >= 50 % of the NCBI tRNA, same amino acid, and the same anticodon
when NCBI gives one (/anticodon seq:), so trnL1/trnL2 and trnS1/trnS2 are distinguished.
False positive = predicted locus overlapping no NCBI tRNA. Input must have been run on bare
sequences (NCBI tRNAs not given to mito_trna.py).

Usage: python3 eval_trna.py <ncbi_gb_dir> <mito_trna_out_root> [label]
  (out root contains one sub-directory per record with trna_report.tsv)
"""
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

GB, ROOT = Path(sys.argv[1]), Path(sys.argv[2])
LABEL = sys.argv[3] if len(sys.argv) > 3 else ROOT.name
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H",
       "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V"}


def truth(rec):
    out = []
    for f in rec.features:
        if f.type != "tRNA":
            continue
        m = re.search(r"tRNA-([A-Z][a-z]{2})", f.qualifiers.get("product", [""])[0])
        if not m or m.group(1) not in AA3:
            continue
        ac = re.search(r"seq:([acgtu]{3})", " ".join(f.qualifiers.get("anticodon", [])), re.I)
        out.append((AA3[m.group(1)], ac.group(1).lower().replace("u", "t") if ac else "",
                    int(f.location.start) + 1, int(f.location.end), f.location.strand))
    return out


def match(t, p):
    aa, ac, s, e, st = t
    paa, pac, ps, pe, pst = p
    return pst == st and paa == aa and (not ac or not pac or ac == pac) and min(e, pe) - max(s, ps) + 1 >= 0.5 * (e - s + 1)


def score(T, P, c):
    used = set()
    for t in T:
        hits = [i for i, p in enumerate(P) if match(t, p)]
        if not hits:
            c["missed"] += 1
            continue
        i = min(hits, key=lambda k: abs(P[k][2] - t[2]) + abs(P[k][3] - t[3]))
        used.add(i)
        d = max(abs(P[i][2] - t[2]), abs(P[i][3] - t[3]))
        c["exact" if d == 0 else "<=3nt" if d <= 3 else ">3nt"] += 1
    c["FP"] += sum(1 for i, p in enumerate(P) if i not in used and not any(
        min(t[3], p[3]) - max(t[2], p[2]) > 0 for t in T))
    c["truth"] += len(T)


res = defaultdict(Counter)
n_rec = 0
for d in sorted(ROOT.iterdir()):
    rep = d / "trna_report.tsv"
    gb = GB / f"{d.name}.gb"
    if not rep.exists() or not gb.exists():
        continue
    n_rec += 1
    T = truth(SeqIO.read(str(gb), "genbank"))
    rows = list(csv.DictReader(rep.open(), delimiter="\t"))
    by_src = defaultdict(list)
    for r in rows:
        aa = {"trnL1": "L", "trnL2": "L", "trnS1": "S", "trnS2": "S"}.get(r["tRNA"], r["tRNA"][3:4])
        st = 1 if r["strand"] == "+" else -1
        for item in r["coords_by_source"].split("; "):
            src, co = item.split(":")
            s, e = map(int, co.split(".."))
            by_src[src].append((aa, r["anticodon"], s, e, st))
        for k in (1, 2, 3):
            if int(r["n_families"]) >= k and not r["conflict"]:
                by_src[f"selected >={k} families"].append((aa, r["anticodon"], int(r["start"]), int(r["end"]), st))
        if r.get("selected", "") == "yes":                       # mito_trna.py default rule
            by_src["selected (default)"].append((aa, r["anticodon"], int(r["start"]), int(r["end"]), st))
    for src in ("MITOS2", "ARWEN", "ARAGORN", "tRNAscan-SE", "selected >=1 families", "selected >=2 families",
                "selected >=3 families", "selected (default)"):
        score(T, by_src.get(src, []), res[src])
print(f"{LABEL}: {n_rec} records, {res['MITOS2']['truth']} NCBI tRNAs")
for src, c in res.items():
    found = c["exact"] + c["<=3nt"] + c[">3nt"]
    print(f"  {src:22s} found {found:4d} ({found / max(c['truth'], 1):5.1%})  exact {c['exact']:4d}  <=3nt {c['<=3nt']:4d}"
          f"  >3nt {c['>3nt']:3d}  missed {c['missed']:4d}  FP {c['FP']:4d}")
