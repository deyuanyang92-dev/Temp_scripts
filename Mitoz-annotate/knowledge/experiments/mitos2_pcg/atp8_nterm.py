"""N-terminal residues and length of the annotated atp8 protein in RefSeq GenBank records (knowledge experiment, 2026-10-20).

Usage: python3 atp8_nterm.py rec1.gb rec2.gb ... > atp8_nterm.tsv   (columns: acc, first 8 aa, length aa)
"""
import sys

from Bio import SeqIO

from pcg_names import pcg_name

print("acc\tfirst8\tlength_aa")
for path in sys.argv[1:]:
    rec = next(SeqIO.parse(path, "genbank"))
    for ft in rec.features:
        if ft.type != "CDS":
            continue
        if pcg_name(ft) == "atp8":
            tr = ft.qualifiers.get("translation", [""])[0]
            print(f"{rec.id}\t{tr[:8]}\t{len(tr)}")
            break
