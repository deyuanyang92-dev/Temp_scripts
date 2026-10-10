"""Compare intron calls: NCBI RefSeq annotation vs real MFannot vs mito_intron.py (full / loo library)."""
import csv
import re
import sys
from collections import Counter
from pathlib import Path

from Bio import SeqIO

B = Path(sys.argv[1])
MARK = re.compile(r"^;\s+G-(\S+?)-I(\d+)\s+(==>|<==)\s+(start|end)\b")


def ncbi_introns(rid):
    """NCBI truth: intron features, plus introns implied by joined CDS (gaps >= 142 nt, no origin wrap)
    when no intron feature overlaps them."""
    r = SeqIO.read(str(B / "fasta" / f"{rid}.gb"), "genbank")
    feats = [(int(f.location.start) + 1, int(f.location.end)) for f in r.features if f.type == "intron"]
    implied = []
    for f in r.features:
        if f.type != "CDS" or len(f.location.parts) < 2:
            continue
        parts = sorted((int(p.start) + 1, int(p.end)) for p in f.location.parts)
        if parts[0][0] == 1 and parts[-1][1] == len(r.seq):
            continue                                    # wraps the origin
        for (a, b), (c, d) in zip(parts, parts[1:]):
            if c - b - 1 >= 142 and not any(min(x[1], c - 1) - max(x[0], b + 1) > 0 for x in feats):
                implied.append((b + 1, c - 1))
    return sorted(feats + implied)


def mfannot_introns(rid):
    mf = B / "mfannot" / f"{rid}.mf"
    if not mf.exists():
        return None
    pending, seen, out, last_end = [], {}, [], 0
    for line in mf.read_text(errors="ignore").splitlines():
        m = MARK.match(line)
        if m:
            pending.append(f"{m.group(1)}-I{m.group(2)}")
            continue
        s = re.match(r"^\s*(\d+)\s+([A-Za-z!]+)", line)
        if s and not line.startswith(";"):
            pos = int(s.group(1))
            for key in pending:
                seen.setdefault(key, []).append(pos)
            pending = []
            last_end = pos + len(s.group(2).replace("!", "")) - 1
    for key in pending:
        seen.setdefault(key, []).append(last_end + 1)
    for key, ps in seen.items():
        if len(ps) >= 2:
            out.append((min(ps[0], ps[1] - 1), max(ps[0], ps[1] - 1)))
    return sorted(out)


def ours_introns(mode, rid):
    rep = B / f"ours_{mode}" / rid / "intron_report.tsv"
    if not rep.exists():
        return None
    return sorted((int(r["intron_start"]), int(r["intron_end"]))
                  for r in csv.DictReader(rep.open(), delimiter="\t") if r.get("intron_start") and r.get("written", "yes") != "no")


def compare(truth, pred):
    """Classify each truth intron against predictions; count unmatched predictions."""
    res, used = Counter(), set()
    for t in truth:
        tl = t[1] - t[0] + 1
        best = None
        for i, p in enumerate(pred):
            ov = max(0, min(t[1], p[1]) - max(t[0], p[0]) + 1)
            if ov >= 0.5 * tl and (best is None or ov > best[1]):
                best = (i, ov)
        if best is None:
            res["missed"] += 1
            continue
        used.add(best[0])
        p = pred[best[0]]
        d = max(abs(p[0] - t[0]), abs(p[1] - t[1]))
        res["exact" if d == 0 else "<=6nt" if d <= 6 else ">6nt"] += 1
    res["extra"] += len(pred) - len(used)
    return res


recs = [line.rstrip("\n").split("\t") for line in (B / "records.tsv").open()]
tools = {"MFannot": mfannot_introns, "ours_full": lambda r: ours_introns("full", r),
         "ours_loo": lambda r: ours_introns("loo", r)}
summary = {k: {"positive": Counter(), "negative": Counter()} for k in tools}
summary["ours_full vs MFannot"] = {"positive": Counter(), "negative": Counter()}
rows = []
for rid, set_, code, org in recs:
    truth = ncbi_introns(rid)
    calls = {k: f(rid) for k, f in tools.items()}
    for k, pred in calls.items():
        if pred is None:
            summary[k][set_]["no output"] += 1
            continue
        c = compare(truth, pred)
        summary[k][set_].update(c)
        summary[k][set_]["truth"] += len(truth)
    if calls["MFannot"] is not None and calls["ours_full"] is not None:
        summary["ours_full vs MFannot"][set_].update(compare(calls["MFannot"], calls["ours_full"]))
    rows.append([rid, set_, org, len(truth)] + [";".join(f"{a}..{b}" for a, b in (calls[k] or [])) for k in tools]
                + [";".join(f"{a}..{b}" for a, b in truth)])
with open(B / "compare.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["record", "set", "organism", "n_ncbi"] + list(tools) + ["ncbi"])
    w.writerows(rows)
for k, v in summary.items():
    for set_, c in v.items():
        print(f"{k:22s} {set_:9s} " + " ".join(f"{x}={c[x]}" for x in
              ("truth", "exact", "<=6nt", ">6nt", "missed", "extra", "no output") if c[x]))
