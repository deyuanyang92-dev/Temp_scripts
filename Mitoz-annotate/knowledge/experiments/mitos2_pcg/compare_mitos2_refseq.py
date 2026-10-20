"""Compare MITOS2 protein-gene boundaries with a RefSeq GenBank record (knowledge experiment, 2026-10-20).

Usage: python3 compare_mitos2_refseq.py <refseq.gb> <mitos_outdir> <transl_table> [auxinfo.json] > rows.tsv
Columns: acc, gene, strand, refseq start..end (1-based), mitos start..end (1-based), d5, d3 (nt; + = MITOS2 shorter at 5',
longer at 3'), RefSeq first codon, MITOS2 first codon, learned start probability of both (auxinfo.json 'start'[code*][gene]),
RefSeq last 3 nt, MITOS2 last 3 nt. result.mitos: columns 5/6 are 0-based inclusive start/stop.
"""
import json
import sys

from Bio import SeqIO

from pcg_names import pcg_name

GENES = ["atp8", "atp6", "nad4l", "nad4", "nad6", "cox1", "cox2", "cox3", "cob", "nad1", "nad2", "nad3", "nad5"]


def main(gb, outdir, code, aux="/opt/mitos_refs/refseq89m/auxinfo.json"):
    rec = next(SeqIO.parse(gb, "genbank"))
    seq = rec.seq
    start_p = json.load(open(aux))["start"][f"{code}*"]
    ref = {}
    for ft in rec.features:
        if ft.type == "CDS" and pcg_name(ft):
            ref[pcg_name(ft)] = (int(ft.location.start) + 1, int(ft.location.end), ft.location.strand)
    mit = {}
    for line in open(f"{outdir}/result.mitos"):
        c = line.rstrip("\n").split("\t")
        if c[1] == "gene":
            mit.setdefault(c[2].split("_")[0].split("-")[0], []).append((int(c[4]) + 1, int(c[5]) + 1, int(c[6])))

    def cod(a, which):
        sub = seq[a[0] - 1:a[1]]
        sub = sub.reverse_complement() if a[2] == -1 else sub
        return str(sub[:3] if which == "start" else sub[-3:])

    print("\t".join(["acc", "gene", "strand", "refseq", "mitos", "d5", "d3", "ref_first", "mitos_first", "p_ref_first",
                     "p_mitos_first", "ref_last3", "mitos_last3", "n_mitos_copies"]))
    for g in GENES:
        a, b = ref.get(g), mit.get(g)
        if not a or not b:
            print("\t".join([rec.id, g, "", str(a), str(b), "", "", "", "", "", "", "", "", ""]))
            continue
        copies = len(b)
        b = max(b, key=lambda x: x[1] - x[0])
        if not 0 < b[1] - b[0] < 10000:                 # coordinate wrapped over the origin
            print("\t".join([rec.id, g, "", str(a), str(b), "wrap", "", "", "", "", "", "", "", str(copies)]))
            continue
        d5 = b[0] - a[0] if a[2] == 1 else a[1] - b[1]
        d3 = b[1] - a[1] if a[2] == 1 else a[0] - b[0]
        rf, mf = cod(a, "start"), cod(b, "start")
        p = lambda c: round(start_p.get(g, {}).get(c, 0.0), 3)
        print("\t".join(map(str, [rec.id, g, a[2], f"{a[0]}..{a[1]}", f"{b[0]}..{b[1]}", d5, d3, rf, mf, p(rf), p(mf),
                                  cod(a, "stop"), cod(b, "stop"), copies])))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), *sys.argv[4:5])
