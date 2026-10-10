"""Robustness: odd records must not crash a run, and good records in the same file must still be annotated.

Run: python3 -m pytest -q tests/test_robustness.py
"""
import csv
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
SCRIPT = HERE.parent / "mito_intron.py"
BIN = str(Path(shutil.which("exonerate") or "/root/miniconda3/envs/mitointron/bin/exonerate").parent)
pytestmark = pytest.mark.skipif(not (Path(BIN) / "exonerate").exists(), reason="exonerate not available")


def rec(rid, seq):
    r = SeqRecord(Seq(seq), id=rid, name=rid, description="odd record")
    r.annotations.update(molecule_type="DNA", topology="circular", organism="Testus oddus")
    return r


def test_odd_records_do_not_break_the_batch(tmp_path):
    good = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
    good.features = [f for f in good.features if f.type == "source"]          # bare sequence
    import random
    random.seed(0)
    odd = [rec("EMPTY1", ""), rec("ALLN1", "N" * 5000), rec("SHORT1", "ACGT" * 12),
           rec("IUPAC1", "".join(random.choice("ACGTRYKMSWN") for _ in range(4000))),
           rec("LOWER1", str(good.seq).lower()[:6000])]
    (tmp_path / "in").mkdir()
    SeqIO.write(odd[:3] + [good] + odd[3:], str(tmp_path / "in" / "mixed.gbf"), "genbank")
    p = subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(tmp_path / "in"), "-o",
                        str(tmp_path / "out"), "--genetic_code", "4", "--bin", BIN, "--timeout", "120"],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-2000:]
    out = list(SeqIO.parse(str(tmp_path / "out" / "mixed.gbf"), "genbank"))
    assert [r.id for r in out] == [r.id for r in odd[:3] + [good] + odd[3:]]       # nothing lost or reordered
    rows = list(csv.DictReader((tmp_path / "out" / "intron_report.tsv").open(), delimiter="\t"))
    written = sorted((int(r["intron_start"]), int(r["intron_end"])) for r in rows
                     if r.get("written") == "yes" and r["record"] == good.id)
    assert written == [(4156, 5008), (7744, 9424)]                                   # good record still exact
