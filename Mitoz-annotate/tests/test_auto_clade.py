"""--clade auto: clade, genetic code and phylum from the GenBank input (batch_mitoz.auto_clade).

Run: python3 -m pytest -q tests/test_auto_clade.py
"""
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import batch_mitoz as bm  # noqa: E402


def gb(path, phylum_lineage, code=None, n=1, empty_taxonomy=False):
    recs = []
    for i in range(n):
        r = SeqRecord(Seq("ATG" * 200), id=f"X{i}.1", name=f"X{i}", description="test")
        r.annotations.update(molecule_type="DNA", topology="circular", organism="Genus species",
                             taxonomy=[] if empty_taxonomy else ["Eukaryota", "Metazoa", *phylum_lineage])
        if code:
            r.features.append(SeqFeature(FeatureLocation(0, 99, 1), type="CDS", qualifiers={"transl_table": [str(code)]}))
        recs.append(r)
    SeqIO.write(recs, str(path), "genbank")
    return path


def test_mitoz_phylum(tmp_path):
    f = gb(tmp_path / "a.gb", ["Spiralia", "Lophotrochozoa", "Annelida", "Polychaeta"], code=5)
    assert bm.auto_clade(f) == ("Annelida-segmented-worms", 5, "Annelida", "")


def test_extended_clade_and_code_from_cds(tmp_path):
    f = gb(tmp_path / "c.gb", ["Cnidaria", "Anthozoa"], code=4)
    assert bm.auto_clade(f)[:3] == ("Cnidaria", 4, "Cnidaria")


def test_other_phylum_uses_animal_db_and_phylum_default_code(tmp_path):
    """No CDS transl_table (e.g. an unannotated GenBank): the phylum's usual code is used."""
    f = gb(tmp_path / "t.gb", ["Ctenophora", "Tentaculata"])
    assert bm.auto_clade(f)[:3] == ("Animal", 4, "Ctenophora")


def test_directory_of_one_phylum(tmp_path):
    gb(tmp_path / "a.gb", ["Ecdysozoa", "Nematoda", "Chromadorea"], code=5, n=2)
    gb(tmp_path / "b.gbk", ["Ecdysozoa", "Nematoda", "Enoplea"], code=5)
    assert bm.auto_clade(tmp_path)[:3] == ("Nematoda", 5, "Nematoda")


def test_mixed_phyla_is_an_error(tmp_path):
    gb(tmp_path / "a.gb", ["Annelida"], code=5)
    gb(tmp_path / "b.gb", ["Chordata", "Vertebrata"], code=2)
    clade, code, phylum, err = bm.auto_clade(tmp_path)
    assert clade is None and "Annelida" in err and "Chordata" in err and "--clade" in err


def test_no_taxonomy_falls_back_to_previous_default(tmp_path):
    """check_mitoz.py 2026-10-21: input with NO taxonomy at all (FASTA, GenBank without a lineage) ran before
    with the old default clade; it still runs, with the old default and a warning, instead of stopping."""
    f = tmp_path / "x.fasta"
    f.write_text(">x\nACGT\n")
    clade, code, phylum, err = bm.auto_clade(f)
    assert (clade, code, phylum) == ("Annelida-segmented-worms", 5, None) and err.startswith("WARNING")
    g = gb(tmp_path / "n.gb", [], empty_taxonomy=True)
    assert bm.auto_clade(g)[:3] == ("Annelida-segmented-worms", 5, None)


def test_records_without_taxonomy_do_not_block_the_known_phylum(tmp_path):
    gb(tmp_path / "a.gb", ["Mollusca", "Gastropoda"], code=5)
    gb(tmp_path / "b.gb", [], empty_taxonomy=True)
    assert bm.auto_clade(tmp_path)[:3] == ("Mollusca", 5, "Mollusca")


# --- review B finding 2: a lineage that names no listed phylum must not fall back to the annelid database ---
def test_nonempty_taxonomy_without_known_phylum_is_an_error(tmp_path):
    """B2: lineage [Eukaryota, Metazoa] only (hand-edited header) -> ERROR, not an annelid MitoZ run."""
    clade, code, phylum, err = bm.auto_clade(gb(tmp_path / "u.gb", [], code=2))
    assert clade is None and code is None and not err.startswith("WARNING")
    assert "--clade" in err and "--genetic_code" in err


def test_nonempty_unmatched_taxonomy_next_to_known_phylum_is_an_error(tmp_path):
    gb(tmp_path / "a.gb", ["Mollusca", "Gastropoda"], code=5)
    gb(tmp_path / "b.gb", ["Strange group"], code=5)
    clade, _, _, err = bm.auto_clade(tmp_path)
    assert clade is None and "Strange group" in err


@pytest.mark.parametrize("token,phylum", [
    ("Vertebrata", "Chordata"), ("Craniata", "Chordata"), ("Tunicata", "Chordata"), ("Cephalochordata", "Chordata"),
    ("Hexapoda", "Arthropoda"), ("Insecta", "Arthropoda"), ("Crustacea", "Arthropoda"),
    ("Chelicerata", "Arthropoda"), ("Myriapoda", "Arthropoda"),
    ("Gastropoda", "Mollusca"), ("Bivalvia", "Mollusca"), ("Cephalopoda", "Mollusca"),
    ("Polyplacophora", "Mollusca"), ("Scaphopoda", "Mollusca"),
    ("Polychaeta", "Annelida"), ("Clitellata", "Annelida"), ("Sipuncula", "Annelida"),
    ("Anthozoa", "Cnidaria"), ("Hydrozoa", "Cnidaria"), ("Scyphozoa", "Cnidaria"), ("Cubozoa", "Cnidaria"),
    ("Demospongiae", "Porifera"), ("Hexactinellida", "Porifera"), ("Calcarea", "Porifera"),
    ("Homoscleromorpha", "Porifera"),
    ("Turbellaria", "Platyhelminthes"), ("Trematoda", "Platyhelminthes"), ("Cestoda", "Platyhelminthes"),
    ("Monogenea", "Platyhelminthes"),
    ("Chromadorea", "Nematoda"), ("Enoplea", "Nematoda")])
def test_class_level_lineage_maps_to_phylum(tmp_path, token, phylum):
    """B2: Geneious / MitoFish style lineages without the phylum token."""
    f = gb(tmp_path / "s.gb", [token], code=5)
    got = bm.auto_clade(f)
    assert got[2] == phylum and got[0] == bm.PHYLUM_CLADE.get(phylum, "Animal") and got[3] == "", got


def test_phylum_token_wins_over_synonym(tmp_path):
    f = gb(tmp_path / "s.gb", ["Mollusca", "Cephalopoda"], code=5)
    assert bm.auto_clade(f)[2] == "Mollusca"


# --- review B finding 3: codes per record ---
def test_disagreeing_genetic_codes_are_an_error(tmp_path):
    """B3: 2 records of code 5 + 3 records of code 2 used to give code 2 to all five, silently."""
    gb(tmp_path / "a.gb", ["Mollusca"], code=5, n=2)
    gb(tmp_path / "b.gb", ["Mollusca"], code=2, n=3)
    clade, code, phylum, err = bm.auto_clade(tmp_path)
    assert clade is None and "2: 3 records" in err and "5: 2 records" in err and "--genetic_code" in err


def test_disagreeing_codes_fine_with_explicit_genetic_code(tmp_path):
    gb(tmp_path / "a.gb", ["Mollusca"], code=5, n=2)
    gb(tmp_path / "b.gb", ["Mollusca"], code=2, n=3)
    assert bm.auto_clade(tmp_path, genetic_code=5)[:3] == ("Mollusca", 5, "Mollusca")


def test_codes_counted_per_record_not_per_cds(tmp_path):
    """One record with many CDS of code 2 and one record with one CDS of code 5 are two records, 1:1."""
    r = gb(tmp_path / "a.gb", ["Mollusca"], code=2)
    recs = list(SeqIO.parse(str(r), "genbank"))
    for i in range(5):
        recs[0].features.append(SeqFeature(FeatureLocation(10 * i, 10 * i + 9, 1), type="CDS",
                                           qualifiers={"transl_table": ["2"]}))
    SeqIO.write(recs, str(r), "genbank")
    gb(tmp_path / "b.gb", ["Mollusca"], code=5)
    err = bm.auto_clade(tmp_path)[3]
    assert "2: 1 records" in err and "5: 1 records" in err


# --- review B findings 12, 13 ---
def test_gbff_is_a_genbank_file_everywhere(tmp_path):
    """B12: auto_clade read .gbff but input sniffing said fasta."""
    f = gb(tmp_path / "a.gbff", ["Mollusca"], code=5)
    assert bm.auto_clade(tmp_path)[2] == "Mollusca"
    assert bm._sniff_input_type(f) == "genbank" and bm._sniff_input_type(tmp_path) == "genbank"
    assert bm.GB_SUFFIXES == tuple(bm._GB_EXTS) or set(bm.GB_SUFFIXES) == set(bm._GB_EXTS)


def test_truncated_genbank_names_the_file(tmp_path):
    """B13: bare ValueError traceback from main() before; now an error message with the file name."""
    bad = tmp_path / "bad2.gb"
    bad.write_text("LOCUS       X                 100 bp    DNA     circular UNK 01-JAN-1980\nFEATURES\nORIGIN\n")
    clade, _, _, err = bm.auto_clade(bad)
    assert clade is None and "bad2.gb" in err and not err.startswith("WARNING")
