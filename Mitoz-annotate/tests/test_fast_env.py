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
    # MitoZ uses the list as a MINIMUM gene length: int(0.5 x median aa x 3); absent gene 0
    assert lengths["COX1"] == "16" and lengths["ND6"] == "7" and lengths["ATP8"] == "0"


def test_profiles_header_regex_and_warning(tmp_path, capsys):
    """A (mito_env.py:209): headers with letters after NC_ or a version suffix were skipped silently."""
    src = tmp_path / "mitoz_profiles"
    (src / "MT_database").mkdir(parents=True)
    (src / "rRNA_CM").mkdir()
    (src / "MT_database" / "Cnidaria_CDS_protein.fa").write_text(
        ">gi_NC_KU507297_ATP6_Bolocera_sp._229_aa\n" + "M" * 200 + "\n"
        ">gi_NC_011820.2_ATP6_Foo_bar_228_aa\n" + "M" * 100 + "\n"
        ">gi_NC_011820.2_ND4L_Foo_bar_98_aa\n" + "M" * 50 + "\n"
        ">gi_NC_0123_ND4_Foo_bar_400_aa\n" + "M" * 400 + "\n"
        ">something_unparsable\n" + "M" * 10 + "\n")
    d = mito_env.mitoz_profiles(src, "Cnidaria", tmp_path / "out", as_clade="Arthropoda")
    lengths = dict(line.split("\t") for line in (d / "CDS_HMM" / "Arthropoda_CDS_length_list").read_text().splitlines())
    assert lengths["ATP6"] == str(int(0.5 * 150 * 3))          # both ATP6 headers parsed: median of 200, 100
    assert lengths["ND4L"] == str(int(0.5 * 50 * 3)) and lengths["ND4"] == str(int(0.5 * 400 * 3))
    err = capsys.readouterr().err
    assert "1 of 5" in err and "mitoz_profiles" in err           # one header did not match


def test_profiles_missing_database_is_a_clear_error_and_leaves_no_link(tmp_path):
    """A (mito_env.py:198-204): a typo in the clade left a dangling symlink, the next call raised FileExistsError."""
    src = _fake_mitoz_profiles(tmp_path)
    for _ in range(2):
        with pytest.raises(FileNotFoundError, match="Cnidaria_typo"):
            mito_env.mitoz_profiles(src, "Cnidaria_typo", tmp_path / "out", as_clade="Arthropoda")
    assert not list((tmp_path / "out").rglob("*_CDS_protein.fa"))


def test_profiles_dangling_symlink_is_replaced(tmp_path):
    src = _fake_mitoz_profiles(tmp_path)
    link = tmp_path / "out" / "Cnidaria" / "MT_database" / "Arthropoda_CDS_protein.fa"
    link.parent.mkdir(parents=True)
    link.symlink_to(tmp_path / "gone.fa")                          # dangling from an earlier failed run
    mito_env.mitoz_profiles(src, "Cnidaria", tmp_path / "out", as_clade="Arthropoda")
    assert link.read_text().startswith(">gi_NC_1_COX1")


SHIPPED = next(iter(sorted((MITOZ_PY.parent.parent / "lib").glob("python*/site-packages/mitoz/profiles"))), None)


@pytest.mark.skipif(SHIPPED is None or not SHIPPED.is_dir(), reason="MitoZ profiles not available")
def test_profile_minimum_lengths_stay_below_mitoz_shipped_lists(tmp_path):
    """A (mito_env.py:186-216): MitoZ drops a genome from control-region detection when a PCG is shorter than
    the list value, so a computed list must not exceed the lists MitoZ ships (they are 0.6-1.0 of the median x 3)."""
    n = 0
    for shipped in sorted((SHIPPED / "CDS_HMM").glob("*_CDS_length_list")):
        clade = shipped.name[:-len("_CDS_length_list")]
        if not (SHIPPED / "MT_database" / f"{clade}_CDS_protein.fa").exists():
            continue
        d = mito_env.mitoz_profiles(SHIPPED, clade, tmp_path / clade, as_clade="Zed")
        mine = dict(line.split("\t") for line in (d / "CDS_HMM" / "Zed_CDS_length_list").read_text().splitlines())
        for line in shipped.read_text().splitlines():
            gene, value = line.split()
            assert int(mine[gene]) <= int(value), (clade, gene, mine[gene], value)
        n += 1
    assert n >= 5


def test_mitos_mode_is_lowercased_and_validated(capsys):
    """C18 (scripts/mitoz_fast.py:54): YZZ_MITOS_MODE=Fast was treated as 'sensitive' without a word."""
    assert mito_env.mitos_mode("Fast") == "fast" and mito_env.mitos_mode(" FULL ") == "full"
    assert mito_env.mitos_mode("sensitive") == "sensitive" and mito_env.mitos_mode(None) == "sensitive"
    assert mito_env.mitos_mode("") == "sensitive"
    assert capsys.readouterr().err == ""
    assert mito_env.mitos_mode("sensetive") == "sensitive"            # unknown -> warning + sensitive
    err = capsys.readouterr().err
    assert "YZZ_MITOS_MODE" in err and "sensetive" in err and "sensitive" in err


def test_mito_trna_rejects_unknown_mitos_mode(tmp_path):
    """C18: any unknown mode silently behaved as 'sensitive'."""
    import types
    import mito_trna
    with pytest.raises(ValueError, match="mitos_mode"):
        mito_trna.mitos2(tmp_path / "g.fa", 5, {"runmitos": "/x/runmitos"}, tmp_path,
                         types.SimpleNamespace(mitos_mode="sensetive", mitos_ref="r", mitos_refdir="d"))


def test_extended_clade_accepted_and_needs_code(tmp_path):
    assert "Animal" in bm.ALL_CLADES and "Cnidaria" in bm.ALL_CLADES
    r = subprocess.run([sys.executable, str(HERE.parent / "batch_mitoz.py"), "run", "-i", str(tmp_path), "-o", str(tmp_path / "o"),
                        "--clade", "Cnidaria"], capture_output=True, text=True)
    assert r.returncode != 0 and "--genetic_code" in (r.stdout + r.stderr)


def test_core_choice_keeps_mitom_steps():
    """--core only picks the annotation core: hybrid (default; MitoZ genes + MITOS2 tRNAs + MITOS2 rescue) or
    mitoz (MitoZ alone: MiTFi tRNAs, no rescue). Reorientation, round-2 annotation, metadata transfer,
    introns are MitoM's own steps (user 2026-10-21) and stay on with either core."""
    p = bm.build_parser()
    a = p.parse_args(["run", "-i", "x", "-o", "y", "--clade", "Mollusca", "--core", "mitoz"])
    assert bm.resolve_tools(a) == "" and (a.trna, a.rescue) == ("mitoz", "no")
    assert (a.which_gene_first, a.if_reannotation_2, a.preserve_phase1_metadata, a.intron_fix) == ("cox1", "yes", "yes", "auto")
    a = p.parse_args(["run", "-i", "x", "-o", "y", "--clade", "Mollusca"])
    assert bm.resolve_tools(a) == "" and (a.core, a.trna, a.rescue) == ("hybrid", "mitos2", "auto")
    a = p.parse_args(["run", "-i", "x", "-o", "y", "--clade", "Animal", "--genetic_code", "5", "--core", "mitoz"])
    assert "MitoZ" in bm.resolve_tools(a)


@pytest.mark.skipif(not mito_env.find("table2asn"), reason="table2asn not available")
def test_run_ends_with_curated_and_submission_files(tmp_path):
    """User 2026-10-21: output must be submittable to NCBI. run = MitoM annotation (04.final_gb) -> curation with
    NCBI names (05.curated) -> per record .fsa/.tbl (BankIt) + .sqn/.val (table2asn) in 06.ncbi_submission."""
    import argparse
    final = tmp_path / "04.final_gb"
    final.mkdir()
    (final / "KY753833.gbf").write_bytes((HERE / "data" / "KY753833_mitoz.gbf").read_bytes())
    args = argparse.Namespace(curate="yes", genetic_code=5, sbt=None)
    merged_from = bm.run_curate_and_submit(args, final, tmp_path)
    assert merged_from == tmp_path / "05.curated"
    rec = SeqIO.read(str(merged_from / "KY753833.gbf"), "genbank")
    assert {f.qualifiers["gene"][0] for f in rec.features if f.type == "rRNA"} == {"rrnS", "rrnL"}
    sub = tmp_path / "06.ncbi_submission"
    rid = rec.id
    assert all((sub / f"{rid}{ext}").exists() for ext in (".fsa", ".tbl", ".sqn", ".val"))
    head = (sub / "submission_summary.tsv").read_text().splitlines()
    assert head[0].startswith("record\tERROR\tWARNING") and head[1].startswith(rid)


def test_doctor_has_required_biopython_row(capsys, monkeypatch):
    """B1: the doctor said `All required programs found` although the scripts could not import Biopython."""
    import Bio
    mito_env.doctor()
    out = capsys.readouterr().out
    assert any(line.lstrip().startswith("OK") and "Biopython" in line and Bio.__version__ in line
               for line in out.splitlines())
    monkeypatch.setitem(sys.modules, "Bio", None)               # import Bio -> ImportError
    assert mito_env.doctor() == 1
    out = capsys.readouterr().out
    assert any(line.lstrip().startswith("!!") and "Biopython" in line for line in out.splitlines())
    assert "pip install biopython" in out and "All required programs found" not in out


@pytest.mark.parametrize("cpus,threads,tasks", [(4, 4, 1), (32, 8, 4), (2, 2, 1), (8, 8, 1), (None, 8, 1)])
def test_default_threads_follow_the_machine(monkeypatch, cpus, threads, tasks):
    """B22: 16 threads x 16 jobs regardless of the machine; now threads = min(8, CPUs), jobs = CPUs // threads."""
    monkeypatch.setattr(bm.os, "cpu_count", lambda: cpus)
    a = bm.build_parser().parse_args(["run", "-i", "x", "-o", "y"])
    assert (a.threads, a.max_tasks) == (threads, tasks)
    a = bm.build_parser().parse_args(["run", "-i", "x", "-o", "y", "--threads", "3", "--max_tasks", "5"])
    assert (a.threads, a.max_tasks) == (3, 5)


def test_prog_name_follows_the_launcher(monkeypatch):
    monkeypatch.delenv("MITOM_LAUNCHER", raising=False)
    assert bm.build_parser().prog == "batch_mitoz.py"
    monkeypatch.setenv("MITOM_LAUNCHER", "1")
    assert bm.build_parser().prog == "mitom"
