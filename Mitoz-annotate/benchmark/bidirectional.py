"""Bidirectional check of mito_intron.py introns against the NCBI records.

Direction 1 (NCBI -> ours): every NCBI intron (intron feature, or the gap of a joined CDS >= 142 nt
that does not wrap the origin) - found exactly / within 6 nt / >6 nt / missed.
Direction 2 (ours -> NCBI): every intron written by mito_intron.py - what NCBI says at that place:
intron feature / joined CDS gap / an uninterrupted CDS spans it (NCBI contradicts) / no annotation
there / record without any annotation. MFannot coordinates are added when <mfannot_dir>/<acc>.mf exists.

Usage: python3 bidirectional.py <ncbi_gb_dir> <ours_dir> <out.tsv> [mfannot_dir ...]
"""
import csv
import re
import sys
from collections import Counter
from pathlib import Path

from Bio import SeqIO

GB, OURS, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
MF_DIRS = [Path(x) for x in sys.argv[4:]]
MARK = re.compile(r"^;\s+G-(\S+?)-I(\d+)\s+(==>|<==)\s+(start|end)\b")


def ncbi_truth(rec):
    L = len(rec)
    out = [(int(f.location.start) + 1, int(f.location.end), "intron feature") for f in rec.features if f.type == "intron"]
    for f in rec.features:
        if f.type != "CDS" or len(f.location.parts) < 2:
            continue
        p = sorted((int(x.start) + 1, int(x.end)) for x in f.location.parts)
        if p[0][0] == 1 and p[-1][1] == L:
            continue
        for a, b in zip(p, p[1:]):
            s, e = a[1] + 1, b[0] - 1
            if e - s + 1 >= 142 and not any(t[0] <= e and t[1] >= s for t in out):
                out.append((s, e, "CDS join gap"))
    return out


def mfannot(acc):
    for d in MF_DIRS:
        mf = d / f"{acc}.mf"
        if not mf.exists():
            continue
        pend, seen, last = [], {}, 0
        for line in mf.read_text(errors="ignore").splitlines():
            m = MARK.match(line)
            if m:
                pend.append(f"{m.group(1)}-I{m.group(2)}")
                continue
            s = re.match(r"^\s*(\d+)\s+([A-Za-z!]+)", line)
            if s and not line.startswith(";"):
                for k in pend:
                    seen.setdefault(k, []).append(int(s.group(1)))
                pend, last = [], int(s.group(1)) + len(s.group(2)) - 1
        for k in pend:
            seen.setdefault(k, []).append(last + 1)
        return [(p[0], p[1] - 1) for p in seen.values() if len(p) >= 2]
    return None


def ov(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


rows, d1, d2 = [], Counter(), Counter()
for f in sorted(GB.glob("*.gb")):
    rec = SeqIO.read(str(f), "genbank")
    rep = OURS / rec.id / "intron_report.tsv"
    if not rep.exists():
        continue
    rr = list(csv.DictReader(rep.open(), delimiter="\t"))
    ours = [(int(r["intron_start"]), int(r["intron_end"]), r["gene"], r["confidence"])
            for r in rr if r.get("intron_start") and r.get("written") == "yes"]
    low = [(int(r["intron_start"]), int(r["intron_end"]), r["gene"], r.get("evidence", ""))
           for r in rr if r.get("intron_start") and r.get("written") == "no"]
    truth = ncbi_truth(rec)
    mf = mfannot(rec.id)
    mf_s = ";".join(f"{a}..{b}" for a, b in mf) if mf is not None else "not run"
    org = rec.annotations.get("organism", "")
    for t in truth:
        hit = max(ours, key=lambda o: ov(t, o), default=None)
        if hit and ov(t, hit) >= 0.5 * (t[1] - t[0] + 1):
            d = max(abs(hit[0] - t[0]), abs(hit[1] - t[1]))
            cls = "exact" if d == 0 else "<=6 nt" if d <= 6 else ">6 nt"
            ours_s = f"{hit[0]}..{hit[1]} ({hit[3]})"
        else:
            lo = [x for x in low if ov(t, x) > 0]
            cls = "missed" + (" (found, not written: " + lo[0][3][:60] + ")" if lo else "")
            ours_s = ""
        d1[cls.split(" (")[0]] += 1
        rows.append(["NCBI->ours", rec.id, org, f"{t[0]}..{t[1]}", t[2], ours_s, cls, mf_s])
    annotated = [x for x in rec.features if x.type not in ("source",)]
    for o in ours:
        if any(ov(o, t) > 0 for t in truth):
            st = next(t[2] for t in truth if ov(o, t) > 0)
        elif not annotated:
            st = "record has no annotation"
        elif any(x.type == "CDS" and len(x.location.parts) == 1 and int(x.location.start) < o[0] and int(x.location.end) >= o[1]
                 for x in rec.features):
            st = "NCBI CDS spans it uninterrupted (contradiction)"
        else:
            st = "no NCBI annotation at this place"
        d2[st] += 1
        rows.append(["ours->NCBI", rec.id, org, f"{o[0]}..{o[1]} {o[2]} ({o[3]})", st, "", "", mf_s])
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["direction", "record", "organism", "intron", "ncbi_status", "ours", "class", "mfannot"])
    w.writerows(rows)
print("NCBI -> ours:", dict(d1))
print("ours -> NCBI:", dict(d2))
