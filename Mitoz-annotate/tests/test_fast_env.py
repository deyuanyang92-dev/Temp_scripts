"""mito_env.py (one place that finds programs), batch_mitoz fast MitoZ launcher, clade check, mitoz_fast.py."""
import os
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import batch_mitoz as bm  # noqa: E402
import mito_env  # noqa: E402

MITOZ_PY = Path(os.environ.get("MITOZ_BIN", "/root/miniconda3/envs/mitoz3.6/bin")) / "python"
MITOS_BIN = mito_env.mitos_bin()
MITOS_REF = mito_env.mitos_refdir()


def fake_exe(d: Path, name: str) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    p = d / name
    p.write_text("#!/bin/sh\nexit 0\n")
    p.chmod(0o755)
    return p


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    """Isolated config file; caches cleared before and after."""
    path = tmp_path / "mito.json"
    monkeypatch.setenv("YZZ_MITO_CONFIG", str(path))
    mito_env.config.cache_clear()
    yield path
    mito_env.config.cache_clear()


def test_config_tools_bin_used(cfg, tmp_path):
    fake = fake_exe(tmp_path / "tools", "exonerate-yzz-test")
    cfg.write_text('{"tools_bin": "%s"}' % fake.parent)
    assert mito_env.find("exonerate-yzz-test") == str(fake)


def test_explicit_wrong_path_is_not_replaced(cfg, tmp_path):
    # 2026-10-20: an explicit --mitos_bin that does not contain runmitos must fail, not silently use another MITOS2
    fake_exe(tmp_path / "tools", "exonerate-yzz-test")
    cfg.write_text('{"tools_bin": "%s"}' % (tmp_path / "tools"))
    assert mito_env.find("exonerate-yzz-test", str(tmp_path / "empty")) is None
    assert mito_env.bin_with("exonerate-yzz-test", explicit=str(tmp_path / "empty")) is None


def test_bin_with_needs_all_programs(cfg, tmp_path):
    d = tmp_path / "half"
    fake_exe(d, "prog-a-yzz")
    assert mito_env.bin_with("prog-a-yzz", "prog-b-yzz") is None
    fake_exe(d, "prog-b-yzz")
    cfg.write_text('{"tools_bin": "%s"}' % d)
    mito_env.config.cache_clear()
    assert mito_env.bin_with("prog-a-yzz", "prog-b-yzz") == str(d)


def test_fast_cmd_conda_and_path():
    fast = "/x/scripts/mitoz_fast.py"
    conda = ["/c/bin/conda", "run", "-n", "mitoz3.6", "--no-capture-output", "mitoz", "annotate", "--clade", "Mollusca"]
    assert bm.fast_mitoz_cmd(conda, fast) == conda[:5] + ["python", fast] + conda[6:]
    assert bm.fast_mitoz_cmd(conda, "") == conda                             # --trna mitoz --circos yes
    assert bm.fast_mitoz_cmd(conda[:6] + ["--version"], fast) == conda[:6] + ["--version"]
    if MITOZ_PY.exists():
        path = [str(MITOZ_PY.parent / "mitoz"), "annotate", "--outprefix", "S1"]
        assert bm.fast_mitoz_cmd(path, fast) == [str(MITOZ_PY), fast, "annotate", "--outprefix", "S1"]


def test_wrong_clade_rejected_before_any_job(tmp_path):
    # 2026-10-20: an unknown --clade only failed per sample after setup (Cnidaria is now an extended clade)
    r = subprocess.run([sys.executable, str(HERE.parent / "batch_mitoz.py"), "run", "-i", str(tmp_path),
                        "-o", str(tmp_path / "out"), "--clade", "Cnidarians"], capture_output=True, text=True)
    assert r.returncode == 2 and "not a supported clade" in r.stderr and not (tmp_path / "out").exists()


@pytest.mark.skipif(not (MITOZ_PY.exists() and MITOS_BIN and MITOS_REF), reason="MitoZ / MITOS2 not available")
def test_mitoz_fast_platynereis_trnas_match_ncbi(tmp_path):
    """MitoZ with MITOS2 tRNAs: Platynereis dumerilii NC_000931, all 22 tRNAs at the NCBI coordinates
    (MiTFi inside MitoZ: 20/22, tRNA-Arg and tRNA-Ala off by 1-2 nt)."""
    ncbi = SeqIO.read(str(HERE / "data" / "NC_000931.gb"), "genbank")
    fa = tmp_path / "S1.fa"
    fa.write_text(f">S1 topology=circular\n{ncbi.seq}\n")
    env = dict(os.environ, YZZ_MITOS_BIN=MITOS_BIN, YZZ_MITOS_REFDIR=MITOS_REF,
               PATH=f"{MITOZ_PY.parent}:{os.environ['PATH']}")
    r = subprocess.run([str(MITOZ_PY), str(HERE.parent / "scripts" / "mitoz_fast.py"), "annotate", "--workdir",
                        str(tmp_path / "w"), "--outprefix", "S1", "--thread_number", "4", "--clade",
                        "Annelida-segmented-worms", "--genetic_code", "5", "--fastafiles", str(fa)],
                       capture_output=True, text=True, env=env, timeout=900)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "MITOS2 fast: 22 tRNA" in r.stdout + r.stderr
    gbf = next((tmp_path / "w").glob("*.result/*.gbf"))
    got = {(int(f.location.start), int(f.location.end), f.location.strand)
           for f in SeqIO.read(str(gbf), "genbank").features if f.type == "tRNA"}
    ref = {(int(f.location.start), int(f.location.end), f.location.strand) for f in ncbi.features if f.type == "tRNA"}
    assert got == ref
    assert not list((tmp_path / "w").glob("*.result/circos.png"))           # circos skipped


def test_wrapped_trna_added_to_final_gbf(tmp_path):
    """2026-10-20 benchmark: trnY across the origin (cox1-first rotation) must reach the final GenBank."""
    import random
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    rnd = random.Random(3)
    r = SeqRecord(Seq("".join(rnd.choice("ACGT") for _ in range(1000))), id="S1", description="")
    r.annotations.update(molecule_type="DNA", topology="circular")
    trna = str((r.seq[979:] + r.seq[:40]).reverse_complement())
    (tmp_path / "final").mkdir()
    SeqIO.write(r, str(tmp_path / "final" / "S1.gbf"), "genbank")
    (tmp_path / "trna_wrap.tsv").write_text(f"S1\tY\tgta\t{trna}\n")
    bm.add_wrapped_trnas(tmp_path / "final", tmp_path / "trna_wrap.tsv")
    out = SeqIO.read(str(tmp_path / "final" / "S1.gbf"), "genbank")
    t = [f for f in out.features if f.type == "tRNA"]
    assert len(t) == 1 and str(t[0].extract(out.seq)) == trna


def _fake_mitoz_profiles(tmp_path):
    src = tmp_path / "mitoz_profiles"
    (src / "MT_database").mkdir(parents=True)
    (src / "rRNA_CM").mkdir()
    (src / "MT_database" / "Cnidaria_CDS_protein.fa").write_text(
        ">gi_NC_1_COX1_Aa_bb_3_aa\nMAAAAAAAAA\n>gi_NC_2_COX1_Cc_dd_3_aa\nMAAAAAAAAAAA\n>gi_NC_2_ND6_Cc_dd_3_aa\nMAAAA\n")
    return src


def test_extended_clade_profiles(tmp_path):
    """2026-10-20: MitoZ ships Cnidaria / flatworm / arrow-worm / all-animal protein databases its --clade
    list does not offer; a profiles dir maps them to the name MitoZ is called with (dev benchmark: 48/48
    records of MitoZ-unsupported phyla annotated, PCG both ends exact 54 % vs MITOS2 49 %)."""
    d = mito_env.mitoz_profiles(_fake_mitoz_profiles(tmp_path), "Cnidaria", tmp_path / "out", as_clade="Arthropoda")
    assert (d / "MT_database" / "Arthropoda_CDS_protein.fa").read_text().startswith(">gi_NC_1_COX1")
    assert (d / "rRNA_CM").is_dir()
    lengths = dict(line.split("\t") for line in (d / "CDS_HMM" / "Arthropoda_CDS_length_list").read_text().splitlines())
    assert lengths["COX1"] == "33" and lengths["ND6"] == "15" and lengths["ATP8"] == "0"   # median aa x 3; absent 0


def test_extended_clade_accepted_and_needs_code(tmp_path):
    assert "Animal" in bm.ALL_CLADES and "Cnidaria" in bm.ALL_CLADES
    r = subprocess.run([sys.executable, str(HERE.parent / "batch_mitoz.py"), "run", "-i", str(tmp_path), "-o", str(tmp_path / "o"),
                        "--clade", "Cnidaria"], capture_output=True, text=True)
    assert r.returncode != 0 and "--genetic_code" in (r.stdout + r.stderr)


def test_engine_mitoz_is_plain_mitoz():
    """--engine mitoz = the original MitoZ path (MiTFi tRNAs, no MITOS2 rescue); mitom (default) combines."""
    p = bm.build_parser()
    a = p.parse_args(["run", "-i", "x", "-o", "y", "--clade", "Mollusca", "--engine", "mitoz"])
    assert bm.resolve_tools(a) == "" and (a.trna, a.rescue) == ("mitoz", "no")
    a = p.parse_args(["run", "-i", "x", "-o", "y", "--clade", "Mollusca"])
    assert bm.resolve_tools(a) == "" and (a.engine, a.trna, a.rescue) == ("mitom", "mitos2", "auto")
    a = p.parse_args(["run", "-i", "x", "-o", "y", "--clade", "Animal", "--genetic_code", "5", "--engine", "mitoz"])
    assert "MitoZ" in bm.resolve_tools(a)
