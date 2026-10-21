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
    (as an NCBI record that does not start at cox1). tail_gene=True: a 28 nt piece of tRNA (cut) at the end."""
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
        r.features.append(SeqFeature(FeatureLocation(L - 30, L - 2, 1), type="tRNA", qualifiers={"gene": ["trnF"]}))
    return r


def with_edge_feature(r, ftype, gene, start, end, strand=1, seq=None):
    """Add a feature (optionally with a chosen sequence) at start..end of the round-1 record."""
    if seq is not None:
        assert len(seq) == end - start
        r.seq = Seq(str(r.seq)[:start] + seq + str(r.seq)[end:])
    r.features.append(SeqFeature(FeatureLocation(start, end, strand), type=ftype, qualifiers={"gene": [gene]}))
    return r


def rand_nt(n):
    return "".join(rnd.choice("ACGT") for _ in range(n))


def write_pair(tmp_path, r1, rotate=True):
    p1, p2 = tmp_path / "r1.gbf", tmp_path / "reor.gbf"
    SeqIO.write(r1, str(p1), "genbank")
    cox1 = next(f for f in r1.features if f.type == "CDS" and f.qualifiers["gene"] == ["COX1"])
    SeqIO.write(bm.rotate_record(r1, int(cox1.location.start)) if rotate else r1, str(p2), "genbank")
    return p1, p2


def test_clean_rotation_needs_no_reannotation(tmp_path):
    assert bm.reorient_check(*write_pair(tmp_path, genome())) == []


def test_gene_at_original_end_needs_reannotation(tmp_path):
    """A TRUNCATED gene touching the ends of the round-1 sequence may be cut there; after rotation it lies
    inside. (A complete gene at the ends is no reason: B5 review.)"""
    probs = bm.reorient_check(*write_pair(tmp_path, genome(tail_gene=True)))
    assert any("end of the round-1 sequence" in p and "trnF" in p for p in probs)


def edge_probs(tmp_path, r, code=5):
    return [x for x in bm.reorient_check(*write_pair(tmp_path, r), code=code) if "end of the round-1 sequence" in x]


def test_complete_trna_at_position_one_does_not_trigger(tmp_path):
    """B5: NCBI vertebrate records start with trnF at position 1; a complete tRNA at the ends is not cut."""
    r = genome()
    L = len(r)
    with_edge_feature(r, "tRNA", "trnF", 0, 68)                       # record starts with a complete tRNA
    with_edge_feature(r, "tRNA", "trnP", L - 69, L - 1, strand=-1)    # and ends with one
    assert edge_probs(tmp_path, r) == []
    assert bm.reorient_check(*write_pair(tmp_path, r)) == []


def test_truncated_trna_at_position_one_triggers(tmp_path):
    r = with_edge_feature(genome(), "tRNA", "trnF", 0, 28)
    assert any("trnF" in x for x in edge_probs(tmp_path, r))


def test_overlong_trna_at_the_end_triggers(tmp_path):
    r = genome()
    L = len(r)
    with_edge_feature(r, "tRNA", "trnL", L - 120, L - 1)               # > 95 nt: merged with a neighbour
    assert any("trnL" in x for x in edge_probs(tmp_path, r))


def test_complete_cds_at_the_ends_does_not_trigger(tmp_path):
    """Start codon at the left end of a + strand CDS, stop codon at its right end (here the CDS at the end of
    the round-1 sequence): both ends complete -> not cut. Minus strand: the edge side flips."""
    r = genome()
    L = len(r)
    with_edge_feature(r, "CDS", "ND3L", 0, 300, seq="ATG" + rand_nt(294) + "TAA")
    with_edge_feature(r, "CDS", "ND2L", L - 330, L - 30, seq="ATG" + rand_nt(294) + "TAG")
    assert edge_probs(tmp_path, r) == []
    # minus strand at the left end: the CDS 3' end (stop) is the sequence end -> reverse complement of ATG..TAA
    r = genome()
    rc = str(Seq("ATG" + rand_nt(294) + "TAA").reverse_complement())
    with_edge_feature(r, "CDS", "ND3L", 0, 300, strand=-1, seq=rc)
    assert edge_probs(tmp_path, r) == []


def test_cds_length_not_multiple_of_three_at_the_end_triggers(tmp_path):
    r = genome()
    L = len(r)
    with_edge_feature(r, "CDS", "ND3L", L - 302, L - 1, seq="ATG" + rand_nt(294) + "TAAG")        # 301 nt
    assert any("ND3L" in x for x in edge_probs(tmp_path, r))


def test_cds_without_start_or_stop_at_the_edge_triggers(tmp_path):
    r = genome()
    with_edge_feature(r, "CDS", "NDX", 0, 300, seq="CCC" + rand_nt(294) + "TAA")          # no start codon at the edge
    assert any("NDX" in x for x in edge_probs(tmp_path, r))
    r = genome()
    L = len(r)
    with_edge_feature(r, "CDS", "NDY", L - 300, L, seq="ATG" + rand_nt(294) + "CCC")      # no stop at the edge
    assert any("NDY" in x for x in edge_probs(tmp_path, r))
    # the interior end is not judged: a stop-less end that is not at the edge does not matter
    r = genome()
    with_edge_feature(r, "CDS", "NDZ", 0, 300, seq="ATG" + rand_nt(294) + "CCC")
    assert edge_probs(tmp_path, r) == []


def test_t_or_ta_stop_tail_counts_as_complete(tmp_path):
    """Mitochondrial genes end in T / TA completed by the poly-A tail: complete, although length % 3 != 0."""
    r = genome()
    L = len(r)
    with_edge_feature(r, "CDS", "COX3L", L - 302, L - 1, seq="ATG" + rand_nt(297) + "T")
    with_edge_feature(r, "CDS", "ND4X", 0, 302, seq="ATG" + rand_nt(297) + "TA")
    assert edge_probs(tmp_path, r) == []


def test_start_codon_follows_the_genetic_code(tmp_path):
    """ATT is a start codon in the invertebrate mitochondrial code (5) but not in the standard code (1)."""
    r = genome()
    with_edge_feature(r, "CDS", "NDX", 0, 300, seq="ATT" + rand_nt(294) + "TAA")
    assert edge_probs(tmp_path, r, code=5) == []
    assert any("NDX" in x for x in edge_probs(tmp_path, r, code=1))
    assert edge_probs(tmp_path, r) == []                                                    # default code is 5


def test_short_rrna_at_the_edge_triggers_and_full_length_does_not(tmp_path):
    r = genome()
    L = len(r)
    with_edge_feature(r, "rRNA", "s-rRNA", L - 500, L - 1)             # < 600 nt: cut
    with_edge_feature(r, "rRNA", "l-rRNA", 0, 1000)                    # >= 900 nt: complete
    probs = edge_probs(tmp_path, r)
    assert any("s-rRNA" in x for x in probs) and not any("l-rRNA" in x for x in probs)
    r = genome()
    with_edge_feature(r, "rRNA", "l-rRNA", 0, 800)                     # < 900 nt: cut
    with_edge_feature(r, "rRNA", "s-rRNA", len(r) - 700, len(r) - 1)   # >= 600 nt: complete
    probs = edge_probs(tmp_path, r)
    assert any("l-rRNA" in x for x in probs) and not any("s-rRNA" in x for x in probs)


def test_absent_atp8_or_atp6_is_not_a_reason(tmp_path):
    """B4: atp8 / atp6 are absent from whole phyla (nematodes, flatworms, many bivalves, anthozoans)."""
    for gone in ("ATP8", "ATP6"):
        probs = bm.reorient_check(*write_pair(tmp_path, genome(names=[n for n in NAMES if n != gone])))
        assert probs == [], (gone, probs)


def test_standard_gene_list_lacks_atp8_and_atp6():
    assert "ATP8" not in bm.STANDARD_GENES and "ATP6" not in bm.STANDARD_GENES
    assert {"COX1", "COX2", "COX3", "CYTB", "ND1", "s-rRNA", "l-rRNA"} <= set(bm.STANDARD_GENES)


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
