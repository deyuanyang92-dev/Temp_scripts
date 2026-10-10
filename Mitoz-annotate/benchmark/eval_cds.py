"""CDS-level agreement of mito_intron.py models with the NCBI record, and NCBI CDS self-consistency.

For every CDS written by mito_intron.py (note contains 'mito_intron') the NCBI CDS of the same
gene region is located (largest overlap, same strand) and compared: exon structure, 5' end,
3' end, protein. The NCBI CDS itself is checked for problems the NCBI validator would also
flag: internal stop codon, illegal start codon, missing stop (neither a complete stop codon
nor a T/TA incomplete stop).

Usage: python3 eval_cds.py <bench_dir> <ours_subdir>   (bench_dir/fasta/<acc>.gb = NCBI records)
"""
import csv
import sys
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable

B, OURS = Path(sys.argv[1]), sys.argv[2]


def parts(loc):
    return sorted((int(p.start) + 1, int(p.end)) for p in loc.parts)


def cds_problems(f, rec, code):
    table = CodonTable.unambiguous_dna_by_id[code]
    nt = str(f.location.extract(rec.seq)).upper()
    probs = []
    prot = str(f.location.extract(rec.seq)[: len(nt) - len(nt) % 3].translate(table=code))
    if "*" in prot.rstrip("*"):
        probs.append("internal_stop")
    if nt[:3] not in table.start_codons:
        probs.append(f"start_{nt[:3]}")
    tail = len(nt) % 3
    if not (nt[-3:] in table.stop_codons and tail == 0) and not (tail and nt[-tail:] in ("T", "TA")):
        probs.append("no_stop")
    return probs


tally, ncbi_prob, rows = Counter(), Counter(), []
for d in sorted((B / OURS).glob("*/")):
    gbf = list(d.glob("*.gbf"))
    ncbi_f = B / "fasta" / f"{d.name}.gb"
    if not gbf or not ncbi_f.exists():
        continue
    ours, ncbi = SeqIO.read(str(gbf[0]), "genbank"), SeqIO.read(str(ncbi_f), "genbank")
    for f in ours.features:
        if f.type != "CDS" or not any("mito_intron" in n for n in f.qualifiers.get("note", [])):
            continue
        code = int(f.qualifiers.get("transl_table", ["4"])[0])
        op = parts(f.location)
        cands = [g for g in ncbi.features if g.type == "CDS" and g.location.strand == f.location.strand]

        def ov(g):
            return sum(max(0, min(b, d_) - max(a, c) + 1) for a, b in op for c, d_ in parts(g.location))
        best = max(cands, key=ov, default=None)
        if best is None or ov(best) == 0:
            tally["no NCBI CDS"] += 1
            continue
        np_ = parts(best.location)
        same = op == np_
        five = (op[0][0] == np_[0][0]) if f.location.strand == 1 else (op[-1][1] == np_[-1][1])
        three = (op[-1][1] == np_[-1][1]) if f.location.strand == 1 else (op[0][0] == np_[0][0])
        intr_o = [(a[1] + 1, b[0] - 1) for a, b in zip(op, op[1:])]
        intr_n = [(a[1] + 1, b[0] - 1) for a, b in zip(np_, np_[1:])]
        probs = cds_problems(best, ncbi, code)
        for p in probs:
            ncbi_prob[p] += 1
        tally["CDS identical" if same else "CDS differs"] += 1
        tally["5' same"] += five
        tally["3' same"] += three
        tally["introns same"] += intr_o == intr_n
        tally["NCBI CDS has problem"] += bool(probs)
        rows.append([d.name, f.qualifiers.get("gene", [""])[0], str(op), str(np_), five, three,
                     intr_o == intr_n, ";".join(probs)])
n = tally["CDS identical"] + tally["CDS differs"]
print(f"{OURS}: {n} CDS compared; " + ", ".join(f"{k} {v}" for k, v in tally.items()))
print("NCBI CDS self-check problems:", dict(ncbi_prob))
with open(B / f"eval_cds_{OURS}.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["record", "gene", "ours", "ncbi", "5p_same", "3p_same", "introns_same", "ncbi_cds_problems"])
    w.writerows(rows)
