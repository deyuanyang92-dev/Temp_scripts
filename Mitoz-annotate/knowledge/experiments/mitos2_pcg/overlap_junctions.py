"""Overlap of adjacent same-strand PCG pairs (atp8/atp6, nad4l/nad4) in RefSeq GenBank records (knowledge experiment, 2026-10-20).

Usage: python3 overlap_junctions.py rec1.gb rec2.gb ... > overlap_junctions.tsv
Columns: acc, pair, strand, overlap_nt, frame_offset (downstream CDS start minus upstream CDS start, mod 3, + strand only),
overlap_sequence (5'->3' on the coding strand), downstream first codon, upstream last codon. Wrap-around and multi-part CDS are skipped.
"""
import sys

from Bio import SeqIO

from pcg_names import pcg_name

print("\t".join(["acc", "pair", "strand", "overlap_nt", "frame_offset", "overlap_seq", "down_first_codon", "up_last_codon"]))
for path in sys.argv[1:]:
    rec = next(SeqIO.parse(path, "genbank"))
    seq = rec.seq
    cds = {}
    for ft in rec.features:
        if ft.type == "CDS" and len(ft.location.parts) == 1 and pcg_name(ft):
            cds[pcg_name(ft)] = (int(ft.location.start) + 1, int(ft.location.end), ft.location.strand)
    for up, dn in (("atp8", "atp6"), ("nad4l", "nad4")):
        if up not in cds or dn not in cds:
            print("\t".join([rec.id, f"{up}/{dn}", "", "missing", "", "", "", ""]))
            continue
        a, b = cds[up], cds[dn]
        if a[2] != b[2]:
            print("\t".join([rec.id, f"{up}/{dn}", "diff", "", "", "", "", ""]))
            continue
        if a[2] == 1:
            lo, hi = b[0], a[1]
            ov = hi - lo + 1
            off = (b[0] - a[0]) % 3
        else:
            lo, hi = max(a[0], b[0]), min(a[1], b[1])
            ov = hi - lo + 1
            off = ""
        if ov <= 0:
            print("\t".join([rec.id, f"{up}/{dn}", str(a[2]), "0", "", "", "", ""]))
            continue
        if ov > 200:                                  # wrap-around artefact
            print("\t".join([rec.id, f"{up}/{dn}", str(a[2]), "skip_long", "", "", "", ""]))
            continue
        sub = seq[lo - 1:hi]
        sub = sub.reverse_complement() if a[2] == -1 else sub
        print("\t".join([rec.id, f"{up}/{dn}", str(a[2]), str(ov), str(off), str(sub), str(sub[:3]), str(sub[-3:])]))
