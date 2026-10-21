"""--if_reannotation_2 auto: re-annotate after reorientation only when the round-1 annotation may be wrong
(user 2026-10-21: always re-annotating wastes resources). batch_mitoz.reorient_check().

Run: python3 -m pytest -q tests/test_reannotation_check.py
"""
import random
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import batch_mitoz as bm  # noqa: E402

rnd = random.Random(3)
NAMES = ["COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND5", "ND4", "ND4L", "ND6", "CYTB", "ND1", "ND2",
         "s-rRNA", "l-rRNA"]


def genome(offset=3000, names=NAMES, tail_gene=False):
    """Circular genome of 15 genes, 600 nt each with 200 nt spacers, starting `offset` nt into the sequence
    (as an NCBI record that does not start at cox1)."""
    L = 15 * 800 + offset + 1000
    seq = "".join(rnd.choice("ACGT") for _ in range(L))
    r = SeqRecord(Seq(seq), id="S1.1", name="S1", description="test")
    r.annotations.update(molecule_type="DNA", topology="circular")
    for i, n in enumerate(names):
        s = offset + i * 800
        kind = "rRNA" if "rRNA" in n else "CDS"
        r.features += [SeqFeature(FeatureLocation(s, s + 600, 1), type="gene", qualifiers={"gene": [n]}),
                       SeqFeature(FeatureLocation(s, s + 600, 1), type=kind, qualifiers={"gene": [n]})]
    if tail_gene:
        r.features.append(SeqFeature(FeatureLocation(L - 70, L - 2, 1), type="tRNA", qualifiers={"gene": ["trnF"]}))
    return r


def write_pair(tmp_path, r1, rotate=True):
    p1, p2 = tmp_path / "r1.gbf", tmp_path / "reor.gbf"
    SeqIO.write(r1, str(p1), "genbank")
    cox1 = next(f for f in r1.features if f.type == "CDS" and f.qualifiers["gene"] == ["COX1"])
    SeqIO.write(bm.rotate_record(r1, int(cox1.location.start)) if rotate else r1, str(p2), "genbank")
    return p1, p2


def test_clean_rotation_needs_no_reannotation(tmp_path):
    assert bm.reorient_check(*write_pair(tmp_path, genome())) == []


def test_gene_at_original_end_needs_reannotation(tmp_path):
    """A gene touching the ends of the round-1 sequence may be cut there; after rotation it lies inside."""
    probs = bm.reorient_check(*write_pair(tmp_path, genome(tail_gene=True)))
    assert any("end of the round-1 sequence" in p for p in probs)


def test_gene_at_end_without_rotation_is_fine(tmp_path):
    """Nothing moved (record already starts at cox1): round 2 would see the same sequence."""
    r = genome(offset=0)
    r.features.append(SeqFeature(FeatureLocation(len(r) - 70, len(r) - 2, 1), type="tRNA", qualifiers={"gene": ["trnF"]}))
    assert bm.reorient_check(*write_pair(tmp_path, r)) == []


def test_missing_standard_genes_need_reannotation(tmp_path):
    probs = bm.reorient_check(*write_pair(tmp_path, genome(names=NAMES[:12] + NAMES[13:])))
    assert any("ND2" in p and "missing" in p for p in probs)


def test_changed_gene_sequence_after_rotation_is_caught(tmp_path):
    p1, p2 = write_pair(tmp_path, genome())
    r = SeqIO.read(str(p2), "genbank")
    nd4 = next(f for f in r.features if f.type == "CDS" and f.qualifiers["gene"] == ["ND4"])
    r.features.remove(nd4)
    r.features.append(SeqFeature(FeatureLocation(int(nd4.location.start) + 3, int(nd4.location.end) + 3, 1),
                                 type="CDS", qualifiers={"gene": ["ND4"]}))
    SeqIO.write(r, str(p2), "genbank")
    assert any("ND4" in p and "differs" in p for p in bm.reorient_check(p1, p2))


def test_gene_across_new_origin_is_caught(tmp_path):
    p1, p2 = write_pair(tmp_path, genome())
    r = SeqIO.read(str(p2), "genbank")
    L = len(r)
    r.features.append(SeqFeature(FeatureLocation(L - 30, L, 1) + FeatureLocation(0, 40, 1), type="tRNA",
                                 qualifiers={"gene": ["trnX"]}))
    SeqIO.write(r, str(p2), "genbank")
    assert any("origin" in p for p in bm.reorient_check(p1, p2))
