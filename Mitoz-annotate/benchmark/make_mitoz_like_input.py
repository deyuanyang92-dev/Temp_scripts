"""Make MitoZ-like input: drop intron features; split CDS/genes -> exon 1 only."""
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
B = Path(sys.argv[1])
for f in sorted((B / "fasta").glob("*.gb")):
    r = SeqIO.read(str(f), "genbank")
    split = {g for ft in r.features if ft.type == "CDS" and len(ft.location.parts) > 1 for g in ft.qualifiers.get("gene", [])}
    feats = []
    for ft in r.features:
        if ft.type == "intron":
            continue
        if len(ft.location.parts) > 1 or (ft.type == "gene" and set(ft.qualifiers.get("gene", [])) & split):
            p = ft.location.parts[0]
            ft.location = FeatureLocation(p.start, p.end, strand=p.strand)
            ft.qualifiers.pop("translation", None)
        feats.append(ft)
    r.features = feats
    d = B / "sim" / r.id
    d.mkdir(parents=True, exist_ok=True)
    SeqIO.write(r, str(d / f"{r.id}.gbf"), "genbank")
