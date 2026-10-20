"""Tests for mito_intron.py against NCBI records with known group I introns.

Metridium senile NC_000933: COX1 and ND5 each contain a group I intron.
Nematostella NC_008164: ND5 is split (CDS join) but the record has no intron feature.
Platynereis dumerilii NC_000931: no introns (negative control).
Library: MFannot protein collections with the record's own species excluded.

Needs exonerate and cmsearch: on PATH, or MITO_INTRON_BIN=/path/to/env/bin.
Run: python3 -m pytest -q tests/test_mito_intron.py
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
SCRIPT = HERE.parent / "mito_intron.py"
BIN = os.environ.get("MITO_INTRON_BIN") or str(Path(shutil.which("exonerate") or
                                                   "/root/miniconda3/envs/mitointron/bin/exonerate").parent)
pytestmark = pytest.mark.skipif(not (Path(BIN) / "exonerate").exists(), reason="exonerate not available")


def run(inp: Path, out: Path, code: int) -> None:
    """MFannot library with the record's own species excluded (no circular validation)."""
    cmd = [sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(inp), "-o", str(out),
           "--genetic_code", str(code), "--bin", BIN, "--exclude_own_species"]
    subprocess.run(cmd, check=True, capture_output=True)


def ncbi_cds(acc):
    rec = SeqIO.read(str(DATA / f"{acc}.gb"), "genbank")
    return {f.qualifiers["gene"][0]: (str(f.location), f.qualifiers["translation"][0])
            for f in rec.features if f.type == "CDS" and "gene" in f.qualifiers}


def intron_models(path: Path):
    rec = SeqIO.read(str(path), "genbank")
    cds = {f.qualifiers["gene"][0]: (str(f.location), f.qualifiers["translation"][0])
           for f in rec.features if f.type == "CDS" and "mito_intron" in str(f.qualifiers.get("note"))}
    introns = [(int(f.location.start) + 1, int(f.location.end), f.qualifiers["note"][0])
               for f in rec.features if f.type == "intron"]
    return cds, introns


def test_metridium_introns_recovered_from_mitoz_like_input(tmp_path):
    """Introns removed and split genes cut to exon 1 (as MitoZ leaves them) -> NCBI models restored."""
    rec = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
    feats = []
    for f in rec.features:
        if f.type == "intron":
            continue
        if len(f.location.parts) > 1 or (f.type == "gene" and f.qualifiers.get("gene", [""])[0] in ("COX1", "ND5")):
            p = f.location.parts[0]
            f.location = FeatureLocation(p.start, p.end, strand=p.strand)
            f.qualifiers.pop("translation", None)
        feats.append(f)
    rec.features = feats
    (tmp_path / "in").mkdir()
    SeqIO.write(rec, str(tmp_path / "in" / "NC_000933.gbf"), "genbank")
    run(tmp_path / "in", tmp_path / "out", 4)
    cds, introns = intron_models(tmp_path / "out" / "NC_000933.gbf")
    ref = ncbi_cds("NC_000933")
    for gene in ("COX1", "ND5"):
        assert cds[gene] == ref[gene], gene               # same CDS coordinates and protein as NCBI
    assert [(s, e) for s, e, _ in introns] == [(4156, 5008), (7744, 9424)]
    assert all("group I intron" in n for _, _, n in introns)


def test_nematostella_missing_intron_added(tmp_path):
    (tmp_path / "in").mkdir()
    shutil.copy(DATA / "NC_008164.gb", tmp_path / "in" / "NC_008164.gbf")
    run(tmp_path / "in", tmp_path / "out", 4)
    cds, introns = intron_models(tmp_path / "out" / "NC_008164.gbf")
    assert cds["ND5"] == ncbi_cds("NC_008164")["ND5"]
    assert [(s, e) for s, e, _ in introns] == [(718, 2337)]


def test_intronless_record_unchanged(tmp_path):
    (tmp_path / "in").mkdir()
    shutil.copy(DATA / "NC_000931.gb", tmp_path / "in" / "NC_000931.gbf")
    run(tmp_path / "in", tmp_path / "out", 5)
    assert (tmp_path / "out" / "NC_000931.gbf").read_bytes() == (DATA / "NC_000931.gb").read_bytes()
    assert (tmp_path / "out" / "intron_report.tsv").read_text().count("\n") == 1   # header only


def mitoz_like(acc: str, tmp_path: Path) -> Path:
    """NCBI record as MitoZ would leave it: intron features removed, split CDS/genes cut to exon 1."""
    rec = SeqIO.read(str(DATA / f"{acc}.gb"), "genbank")
    split = {g for f in rec.features if f.type == "CDS" and len(f.location.parts) > 1 for g in f.qualifiers.get("gene", [])}
    feats = []
    for f in rec.features:
        if f.type == "intron":
            continue
        if len(f.location.parts) > 1 or (f.type == "gene" and set(f.qualifiers.get("gene", [])) & split):
            p = f.location.parts[0]
            f.location = FeatureLocation(p.start, p.end, strand=p.strand)
            f.qualifiers.pop("translation", None)
        feats.append(f)
    rec.features = feats
    d = tmp_path / "in"
    d.mkdir()
    SeqIO.write(rec, str(d / f"{acc}.gbf"), "genbank")
    return d


def report(out: Path):
    import csv
    return [r for r in csv.DictReader((out / "intron_report.tsv").open(), delimiter="\t") if r["intron_start"]]


def test_mini_exon_plakinastrella(tmp_path):
    """Sponge COX1: two introns separated by a 9-nt exon (NCBI 11936..11944) must not be merged."""
    run(mitoz_like("NC_010217", tmp_path), tmp_path / "out", 4)
    cds, introns = intron_models(tmp_path / "out" / "NC_010217.gbf")
    assert len(introns) == 2
    rec = SeqIO.read(str(tmp_path / "out" / "NC_010217.gbf"), "genbank")
    cox1 = [f for f in rec.features if f.type == "CDS" and "mito_intron" in str(f.qualifiers.get("note"))][0]
    assert (11936, 11944) in [(int(p.start) + 1, int(p.end)) for p in cox1.location.parts]
    assert introns[1][:2] == (11945, 13062)


def test_conflicting_model_not_written(tmp_path):
    """Montipora: a model fusing an ND4-region exon with nad5 exon 2 must not overwrite MitoZ's ND4
    (since 2026-10-21 that model is replaced by a re-aligned one, see test_next_chain_when_best_conflicts)."""
    run(mitoz_like("NC_006902", tmp_path), tmp_path / "out", 4)
    rows = report(tmp_path / "out")
    assert all(r["written"] == "no" for r in rows if "overlaps existing ND4" in r["evidence"])
    out_rec = SeqIO.read(str(tmp_path / "out" / "NC_006902.gbf"), "genbank") \
        if (tmp_path / "out" / "NC_006902.gbf").exists() else None
    if out_rec is not None:
        nd4 = [f for f in out_rec.features if f.type == "CDS" and f.qualifiers.get("gene") == ["ND4"]]
        assert nd4 and str(nd4[0].location) == ncbi_cds("NC_006902")["ND4"][0]   # MitoZ's ND4 untouched


def test_long_nd5_three_prime_end(tmp_path):
    """NC_056771 ND5 is 724 aa: the 3' end must reach the first in-frame stop (NCBI 3902)."""
    run(mitoz_like("NC_056771", tmp_path), tmp_path / "out", 4)
    cds, _ = intron_models(tmp_path / "out" / "NC_056771.gbf")
    assert cds["ND5"][0] == ncbi_cds("NC_056771")["ND5"][0]


def test_submission_qualifiers(tmp_path):
    """Every CDS written has /product; T/TA-terminated CDS carry /transl_except and the polyA note."""
    shutil.copy(DATA / "NC_008164.gb", tmp_path / "x.gbf")
    (tmp_path / "in").mkdir()
    shutil.move(str(tmp_path / "x.gbf"), tmp_path / "in" / "NC_008164.gbf")
    run(tmp_path / "in", tmp_path / "out", 4)
    rec = SeqIO.read(str(tmp_path / "out" / "NC_008164.gbf"), "genbank")
    ours = [f for f in rec.features if f.type == "CDS" and "mito_intron" in str(f.qualifiers.get("note"))]
    assert ours and all(f.qualifiers.get("product") for f in ours)
    for f in ours:
        tail = len(f.location.extract(rec.seq)) % 3
        assert bool(tail) == ("transl_except" in f.qualifiers)


def test_ref_gb_names_and_filters():
    """--ref_gb: CDS translations become homologs named like the MFannot library (nad5, cox1...)."""
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    lib = mi.load_ref_gb([str(DATA / "NC_000933.gb")], 4)
    genes = sorted({p.gene for p in lib.values()})
    assert "cox1" in genes and "nad5" in genes and "nad4L" in genes      # MFannot spelling
    assert all("*" not in p.seq for p in lib.values())
    assert {p.species for p in lib.values()} == {"metridium senile"}


def test_partial_cds_bare_fragment(tmp_path):
    """Travisia COX1 fragment, all annotation removed: the gene runs off the 5' end, so the CDS is
    5'-partial (no start codon) - the intron (NCBI 236..1100, group II) must still be written."""
    from Bio.SeqFeature import BeforePosition
    rec = SeqIO.read(str(DATA / "LC670763.gb"), "genbank")
    rec.features = [f for f in rec.features if f.type == "source"]
    (tmp_path / "in").mkdir()
    SeqIO.write(rec, str(tmp_path / "in" / "LC670763.gbf"), "genbank")
    run(tmp_path / "in", tmp_path / "out", 5)
    _, introns = intron_models(tmp_path / "out" / "LC670763.gbf")
    assert [(s, e) for s, e, _ in introns] == [(236, 1100)] and "group II" in introns[0][2]
    cds = [f for f in SeqIO.read(str(tmp_path / "out" / "LC670763.gbf"), "genbank").features if f.type == "CDS"][0]
    assert isinstance(cds.location.parts[0].start, BeforePosition)


def test_next_chain_when_best_conflicts(tmp_path):
    """2026-10-21 dev benchmark NC_006902.1: the best exonerate chain puts nad5 exon 1 on ND4 (nad4/nad5 are
    homologous) and is refused; the next chain that does not overlap another gene gives NCBI's ND5."""
    run(mitoz_like("NC_006902", tmp_path), tmp_path / "out", 4)
    cds, _ = intron_models(tmp_path / "out" / "NC_006902.gbf")
    assert cds["ND5"][0] == ncbi_cds("NC_006902")["ND5"][0]
