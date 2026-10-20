"""Tests for benchmark/mito_sample.py (Metazoa-wide NC_ benchmark set). No network.

Run: python3 -m pytest -q tests/test_mito_sample.py
"""
import sys
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "benchmark"))
import mito_sample as ms  # noqa: E402


def rec(tables=(), comment="", structured=None, topology="circular", genes=("cox1",), desc="x mitochondrion, complete genome"):
    r = SeqRecord(Seq("A" * 300), id="NC_000001.1", description=desc)
    r.annotations.update({"molecule_type": "DNA", "topology": topology})
    if comment:
        r.annotations["comment"] = comment
    if structured:
        r.annotations["structured_comment"] = structured
    for i, t in enumerate(tables):
        r.features.append(SeqFeature(FeatureLocation(i * 10, i * 10 + 9, 1), type="CDS", qualifiers={"transl_table": [str(t)]}))
    for g in genes:
        r.features.append(SeqFeature(FeatureLocation(100, 200, 1), type="gene", qualifiers={"gene": [g]}))
    return r


def test_majority_code():
    assert ms.majority_code(rec(tables=(5, 5, 5, 4)), "Annelida") == (5, False)
    assert ms.majority_code(rec(), "Annelida") == (5, True)
    assert ms.majority_code(rec(), "Chordata") == (2, True)


def test_evidence_tier():
    assert ms.evidence_tier(rec(structured={"Assembly-Data": {"Sequencing Technology": "PacBio Sequel"}})) == "A-longread"
    assert ms.evidence_tier(rec(structured={"Assembly-Data": {"Sequencing Technology": "Oxford Nanopore MinION"}})) == "A-longread"
    assert ms.evidence_tier(rec(comment="Assembled from RNA-seq reads")) == "A-transcriptome"
    assert ms.evidence_tier(rec(comment="Sanger")) == "B"
    assert ms.evidence_tier(None) == "B"


def test_annotation_method():
    assert ms.annotation_method(rec(), rec(comment="annotated with MITOS WebServer")) == "MITOS"
    assert ms.annotation_method(rec(), rec(comment="Annotated by MitoZ v2.4")) == "MitoZ"
    s = {"Genome-Annotation-Data": {"Annotation Method": "manual curation"}}
    assert ms.annotation_method(rec(), rec(structured=s)) == "other:manual curation"
    assert ms.annotation_method(rec(), None) == "unknown"


def test_structure_flags():
    assert "linear" in ms.structure_flags(rec(topology="linear"))
    assert "no_cox1" in ms.structure_flags(rec(genes=("cob",)))
    assert "multi_chromosome" in ms.structure_flags(rec(desc="Pediculus mitochondrion minichromosome 3, complete sequence"))
    assert ms.structure_flags(rec()) == []


def rows(spec):
    """spec: list of (order, genus, tier, n)."""
    out = []
    for order, genus, tier, n in spec:
        for i in range(n):
            out.append({"accession": f"NC_{order}{genus}{i}", "order": order, "genus": genus, "evidence_tier": tier})
    return out


def test_pick_order_cap_then_fill():
    r = rows([("O1", "G1", "B", 10), ("O2", "G2", "B", 10), ("O3", "G3", "A-longread", 10)])
    got = ms.pick(r, 12, seed=1)
    assert len(got) == 12
    r = rows([(f"O{i}", f"G{i}", "B", 3) for i in range(10)])
    got = ms.pick(r, 12, seed=1)
    assert len(got) == 12
    per = {}
    for x in got:
        per[x["order"]] = per.get(x["order"], 0) + 1
    assert max(per.values()) <= 2


def test_pick_prefers_tier_a():
    r = rows([("O1", "G1", "B", 5), ("O1", "G2", "A-transcriptome", 1)])
    got = ms.pick(r, 2, seed=3)
    assert any(x["evidence_tier"].startswith("A") for x in got)


def test_pick_small_pool():
    r = rows([("O1", "G1", "B", 3)])
    assert len(ms.pick(r, 6, seed=1)) == 3


def test_split_no_genus_leak():
    r = rows([("O1", f"G{i}", "B", 1 + i % 3) for i in range(6)]) + rows([("O2", "G9", "B", 4)])
    for x in r:
        x["phylum"] = "Annelida"
    ms.split_by_genus(r, seed=7)
    dev = {x["genus"] for x in r if x["split"] == "dev"}
    test = {x["genus"] for x in r if x["split"] == "test"}
    assert dev and test and not dev & test
    again = [dict(x, split="") for x in r]
    ms.split_by_genus(again, seed=7)
    assert [x["split"] for x in again] == [x["split"] for x in r]


def test_mitoz_clade():
    assert ms.MITOZ_CLADE["Annelida"] == "Annelida-segmented-worms"
    assert ms.MITOZ_CLADE["Cnidaria"] is None


def test_ranks_from_lineage_xml():
    xml = """<TaxaSet><Taxon><TaxId>6340</TaxId><ScientificName>Platynereis dumerilii</ScientificName><LineageEx>
    <Taxon><TaxId>6340</TaxId><ScientificName>Annelida</ScientificName><Rank>phylum</Rank></Taxon>
    <Taxon><TaxId>1</TaxId><ScientificName>Polychaeta</ScientificName><Rank>class</Rank></Taxon>
    <Taxon><TaxId>2</TaxId><ScientificName>Phyllodocida</ScientificName><Rank>order</Rank></Taxon>
    <Taxon><TaxId>3</TaxId><ScientificName>Nereididae</ScientificName><Rank>family</Rank></Taxon>
    <Taxon><TaxId>4</TaxId><ScientificName>Platynereis</ScientificName><Rank>genus</Rank></Taxon>
    </LineageEx></Taxon></TaxaSet>"""
    got = ms.ranks_from_lineage_xml(xml)
    assert got["6340"] == {"phylum": "Annelida", "class": "Polychaeta", "order": "Phyllodocida",
                           "family": "Nereididae", "genus": "Platynereis"}
