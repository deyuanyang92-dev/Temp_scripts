"""Tests for mito_rrna_hmm.py (last-resort rRNA search with per-phylum profile HMMs).

Run: python3 -m pytest -q tests/test_mito_rrna_hmm.py   (the build/search test needs mafft, hmmbuild, nhmmer)
"""
import random
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_env  # noqa: E402
import mito_rrna_hmm as rh  # noqa: E402

rnd = random.Random(5)


def rand(n):
    return "".join(rnd.choice("ACGT") for _ in range(n))


def rec(seq, feats=(), taxonomy=("Eukaryota", "Metazoa", "Ctenophora", "Tentaculata")):
    r = SeqRecord(Seq(seq), id="S1", description="")
    r.annotations.update(molecule_type="DNA", topology="circular", taxonomy=list(taxonomy))
    for kind, s, e, strand, gene in feats:
        loc = FeatureLocation(s, e, strand)
        r.features += [SeqFeature(loc, type="gene", qualifiers={"gene": [gene]}),
                       SeqFeature(loc, type=kind, qualifiers={"gene": [gene]})]
    return r


def hit(env_from, env_to, hmm_from, hmm_to, strand=1, score=100.0):
    return rh.Hit(env_from, env_to, strand, hmm_from, hmm_to, score)


def test_missing_rrnas_mitoz_and_ncbi_names():
    assert rh.missing(rec(rand(100), [("rRNA", 0, 50, 1, "s-rRNA")])) == ["rrnL"]
    assert rh.missing(rec(rand(100), [("rRNA", 0, 50, 1, "rrnL")])) == ["rrnS"]
    assert rh.missing(rec(rand(100), [("rRNA", 0, 40, 1, "rrnS"), ("rRNA", 50, 90, 1, "l-rRNA")])) == []


def test_phylum_from_taxonomy_or_override():
    lib = {"Ctenophora", "Nematoda"}
    assert rh.phylum_of(rec(rand(10)), lib) == "Ctenophora"
    assert rh.phylum_of(rec(rand(10), taxonomy=["Eukaryota", "Metazoa", "Porifera"]), lib) is None
    assert rh.phylum_of(rec(rand(10)), lib, override="Nematoda") == "Nematoda"


def test_fragments_of_one_gene_merged():
    """Two hits covering consecutive parts of the model (1-400, 420-900) in genome order are one gene."""
    best = rh.best_span([hit(1000, 1400, 1, 400), hit(1450, 1950, 420, 900, score=80)])
    assert (best.start, best.end) == (1000, 1950)


def test_duplicate_copies_not_merged():
    """User 2026-10-21: mitochondrial genes can be duplicated. Two hits over the same model columns are
    two copies, never one gene: only the better copy is kept."""
    best = rh.best_span([hit(1000, 1900, 1, 900, score=90), hit(2000, 2900, 1, 900, score=120)])
    assert (best.start, best.end) == (2000, 2900)


def test_model_order_must_match_genome_order():
    """Plus strand: the later genome piece must cover later model columns; minus strand the reverse."""
    assert (rh.best_span([hit(1000, 1400, 500, 900), hit(1450, 1850, 1, 400, score=80)]).end) == 1400
    m = rh.best_span([hit(1000, 1400, 500, 900, strand=-1), hit(1450, 1850, 1, 400, strand=-1, score=80)])
    assert (m.start, m.end) == (1000, 1850)


def test_far_pieces_not_merged():
    assert rh.best_span([hit(1000, 1400, 1, 400), hit(2000, 2500, 420, 900, score=80)]).end == 1400


def _lib(tmp_path, median=900):
    tmp_path.mkdir(parents=True, exist_ok=True)
    (tmp_path / "manifest.tsv").write_text("phylum\tgene\tn_ref\tn_used\tmedian_len\n"
                                           f"Ctenophora\trrnL\t30\t30\t{median}\n")
    (tmp_path / "Ctenophora.rrnL.hmm").write_text("x")
    return tmp_path


def test_search_adds_rrna_with_note(tmp_path):
    lib = _lib(tmp_path / "lib")
    r = rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")])
    rows = rh.search_record(r, lib, run=lambda hmm, fa: [hit(2000, 2880, 1, 880)])
    assert [(x["gene"], x["status"]) for x in rows] == [("l-rRNA", "added")]
    f = [f for f in r.features if f.type == "rRNA" and f.qualifiers["gene"] == ["l-rRNA"]][0]
    assert (int(f.location.start), int(f.location.end)) == (2000, 2880)
    assert f.qualifiers["product"] == ["16S ribosomal RNA"] and "profile HMM" in f.qualifiers["note"][0]
    assert any(g.type == "gene" and g.qualifiers["gene"] == ["l-rRNA"] for g in r.features)


def test_search_rejects_hit_on_existing_gene_short_or_long(tmp_path):
    lib = _lib(tmp_path / "lib")
    for h, why in ((hit(200, 1300, 1, 880), "existing"), (hit(2000, 2200, 1, 200), "short"),
                   (hit(1600, 3200, 1, 880), "long")):
        r = rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")])
        rows = rh.search_record(r, lib, run=lambda hmm, fa, h=h: [h])
        assert rows[0]["status"].startswith("rejected") and why in rows[0]["status"], (why, rows)
        assert not any(f.qualifiers.get("gene") == ["l-rRNA"] for f in r.features)


def test_search_skips_when_present_or_no_library(tmp_path):
    lib = _lib(tmp_path / "lib")
    never = lambda hmm, fa: (_ for _ in ()).throw(AssertionError("nhmmer must not run"))  # noqa: E731
    full = rec(rand(3000), [("rRNA", 0, 800, 1, "s-rRNA"), ("rRNA", 900, 1900, 1, "l-rRNA")])
    assert rh.search_record(full, lib, run=never) == []
    other = rec(rand(3000), [("rRNA", 0, 800, 1, "s-rRNA")], taxonomy=["Eukaryota", "Metazoa", "Porifera"])
    rows = rh.search_record(other, lib, run=never)
    assert rows[0]["status"] == "skipped: no library for this phylum"


def test_search_dir_in_place(tmp_path):
    lib = _lib(tmp_path / "lib")
    d = tmp_path / "final"
    d.mkdir()
    SeqIO.write(rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")]),
                str(d / "S1.gbf"), "genbank")
    rows = rh.search_dir(d, lib, run=lambda hmm, fa: [hit(2000, 2880, 1, 880)])
    assert [(x["record"], x["gene"], x["status"]) for x in rows] == [("S1", "l-rRNA", "added")]
    genes = [f.qualifiers["gene"][0] for f in SeqIO.read(str(d / "S1.gbf"), "genbank").features if f.type == "rRNA"]
    assert sorted(genes) == ["l-rRNA", "s-rRNA"]


def _mutate(s, rate):
    return "".join(rnd.choice("ACGT".replace(c, "")) if rnd.random() < rate else c for c in s)


@pytest.mark.skipif(not all(mito_env.find(t) for t in ("mafft", "hmmbuild", "nhmmer")),
                    reason="mafft / hmmbuild / nhmmer not found")
def test_build_and_search_real_tools(tmp_path):
    """End to end: library from 6 reference genomes of one phylum, then a new genome whose 16S differs by
    ~15 % from every reference is found by nhmmer."""
    anc = rand(900)
    refs = []
    for i in range(6):
        r = SeqRecord(Seq(rand(400) + _mutate(anc, 0.12) + rand(400)), id=f"R{i}", description="")
        r.annotations.update(molecule_type="DNA", topology="circular", organism=f"Genus{i} sp{i}",
                             taxonomy=["Eukaryota", "Metazoa", "Ctenophora"])
        r.features.append(SeqFeature(FeatureLocation(400, 400 + len(r.seq) - 800, 1), type="rRNA",
                                     qualifiers={"product": ["16S ribosomal RNA"]}))
        refs.append(r)
    gb = tmp_path / "Ctenophora.gb"
    SeqIO.write(refs, str(gb), "genbank")
    lib = tmp_path / "lib"
    rh.build([gb], lib)
    assert (lib / "Ctenophora.rrnL.hmm").exists()
    q = rec(rand(1500) + _mutate(anc, 0.15) + rand(1500), [("CDS", 0, 1400, 1, "COX1")])
    rows = rh.search_record(q, lib)
    added = [x for x in rows if x["gene"] == "l-rRNA"]
    assert added and added[0]["status"] == "added", rows
    f = [f for f in q.features if f.type == "rRNA"][0]
    assert abs(int(f.location.start) - 1500) < 60 and abs(int(f.location.end) - 2400) < 60
