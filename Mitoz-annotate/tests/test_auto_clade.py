"""--clade auto: clade, genetic code and phylum from the GenBank input (batch_mitoz.auto_clade).

Run: python3 -m pytest -q tests/test_auto_clade.py
"""
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import batch_mitoz as bm  # noqa: E402


def gb(path, phylum_lineage, code=None, n=1):
    recs = []
    for i in range(n):
        r = SeqRecord(Seq("ATG" * 200), id=f"X{i}.1", name=f"X{i}", description="test")
        r.annotations.update(molecule_type="DNA", topology="circular", organism="Genus species",
                             taxonomy=["Eukaryota", "Metazoa", *phylum_lineage])
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
    """check_mitoz.py 2026-10-21: GenBank without taxonomy (and FASTA) ran before with the old default
    clade; it still runs, with the old default and a warning, instead of stopping."""
    f = tmp_path / "x.fasta"
    f.write_text(">x\nACGT\n")
    clade, code, phylum, err = bm.auto_clade(f)
    assert (clade, code, phylum) == ("Annelida-segmented-worms", 5, None) and err.startswith("WARNING")
    g = gb(tmp_path / "n.gb", [])
    assert bm.auto_clade(g)[:3] == ("Annelida-segmented-worms", 5, None)


def test_records_without_taxonomy_do_not_block_the_known_phylum(tmp_path):
    gb(tmp_path / "a.gb", ["Mollusca", "Gastropoda"], code=5)
    gb(tmp_path / "b.gb", [])
    assert bm.auto_clade(tmp_path)[:3] == ("Mollusca", 5, "Mollusca")
