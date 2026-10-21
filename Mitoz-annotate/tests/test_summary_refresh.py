"""annotation-summary.tsv counts the final GenBank files, after rescue / rRNA HMM / curation
(2026-10-21 end-to-end run NC_016117.1: summary said PCG 7 rRNA 0, final file had 11 CDS and 2 rRNA).

Run: python3 -m pytest -q tests/test_summary_refresh.py
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


def test_counts_taken_from_final_files(tmp_path):
    (tmp_path / "annotation-summary.tsv").write_text(
        "sample_name\tPCG\ttRNA\trRNA\tcox1_orientation_ok\tbody_source\tmetadata_mode\n"
        "NC_016117.1\t7\t0\t0\tyes\tphase2_mitoz\tmetadata_transferred\n"
        "OTHER\t13\t22\t2\tyes\tphase2_mitoz\tmetadata_transferred\n")
    d = tmp_path / "05.curated"
    d.mkdir()
    r = SeqRecord(Seq("A" * 3000), id="NC_016117.1", description="")
    r.annotations["molecule_type"] = "DNA"
    r.features = [SeqFeature(FeatureLocation(i * 100, i * 100 + 90, 1), type=t, qualifiers={"gene": [f"g{i}"]})
                  for i, t in enumerate(["CDS"] * 11 + ["rRNA"] * 2 + ["tRNA"])]
    SeqIO.write(r, str(d / "NC_016117.1.gbf"), "genbank")
    bm.refresh_summary(tmp_path, d)
    rows = [x.split("\t") for x in (tmp_path / "annotation-summary.tsv").read_text().splitlines()]
    assert rows[1] == ["NC_016117.1", "11", "1", "2", "yes", "phase2_mitoz", "metadata_transferred"]
    assert rows[2][:4] == ["OTHER", "13", "22", "2"]          # no final file: row kept as it was
