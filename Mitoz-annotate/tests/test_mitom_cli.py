"""The `mitom` command: one entry point for every MitoM step.

Run: python3 -m pytest -q tests/test_mitom_cli.py
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_env  # noqa: E402

MITOM = HERE.parent / "mitom"


def run(*a, py=None, env=None):
    return subprocess.run([py or sys.executable, str(MITOM), *a], capture_output=True, text=True, timeout=120,
                          env=env)


def _tools_present():
    return all(mito_env.find(n) for n in mito_env.REQUIRED) and mito_env.mitos_refdir()


def test_help_lists_commands():
    p = run("-h")
    assert p.returncode == 0
    for c in ("run", "doctor", "curate", "intron", "trna", "rescue", "rrna-hmm", "validate"):
        assert c in p.stdout, c


def test_run_help_is_the_workflow():
    p = run("run", "-h")
    assert p.returncode == 0 and "--rrna_hmm" in p.stdout and "--clade" in p.stdout


def test_step_help_forwarded():
    assert "build" in run("rrna-hmm", "-h").stdout
    assert run("curate", "-h").returncode == 0


def test_unknown_command():
    p = run("frobnicate")
    assert p.returncode == 2 and "unknown command" in p.stderr


def test_version():
    p = run("--version")
    assert p.returncode == 0 and p.stdout.startswith("MitoM ")


@pytest.mark.skipif(not _tools_present(), reason="required programs not installed here")
def test_doctor_runs_and_exits_zero():
    """B23 review: test_mitom_cli never ran `mitom doctor`; the doctor now also checks Biopython (B1)."""
    p = run("doctor")
    assert p.returncode == 0, p.stdout[-800:] + p.stderr[-800:]
    assert "Biopython" in p.stdout and "All required programs found" in p.stdout


def test_launcher_sets_prog_name():
    """B22: mitom run -h shows `mitom run`, not batch_mitoz.py."""
    p = run("run", "-h")
    assert "usage: mitom run" in p.stdout


@pytest.fixture
def venv_without_bio(tmp_path):
    """A python that cannot import Biopython (venv without system site packages)."""
    d = tmp_path / "venv"
    subprocess.run([sys.executable, "-m", "venv", "--without-pip", str(d)], check=True, capture_output=True)
    py = d / "bin" / "python"
    assert subprocess.run([str(py), "-c", "import Bio"], capture_output=True).returncode != 0
    return str(py)


def test_launcher_without_biopython_gives_hint(venv_without_bio):
    """B1/B18 review: Biopython missing and no conda env offers it -> one clear message, exit 2 (not a
    traceback or `ERROR: Biopython not installed` from deep inside batch_mitoz.py)."""
    env = dict(os.environ, YZZ_MITO_NO_SCAN="1")
    p = run("run", "-h", py=venv_without_bio, env=env)
    assert p.returncode == 2
    assert f"Biopython not found for {venv_without_bio}" in p.stderr and "install.sh" in p.stderr


def test_launcher_reexecs_with_conda_env_python_that_has_biopython(venv_without_bio, tmp_path):
    """B1/B18 review: the launcher runs under whatever python3 is first on PATH (Miniforge python in Docker,
    an activated conda env); it must re-exec with the python of a conda env that imports Biopython."""
    root = tmp_path / "fakeconda"
    (root / "envs" / "mitotools" / "bin").mkdir(parents=True)
    (root / "envs" / "mitotools" / "bin" / "python").symlink_to(sys.executable)
    env = {k: v for k, v in os.environ.items() if k not in ("CONDA_EXE", "MAMBA_EXE", "CONDA_PREFIX",
                                                              "MAMBA_ROOT_PREFIX", "CONDA_ROOT")}
    env["CONDA_ROOT"] = str(root)
    p = run("run", "-h", py=venv_without_bio, env=env)
    assert p.returncode == 0 and "usage: mitom run" in p.stdout, p.stderr[-600:]


def test_doctor_without_biopython_still_runs_and_flags_it(venv_without_bio):
    """doctor is the command that reports a missing Biopython, so the launcher lets it through."""
    env = dict(os.environ, YZZ_MITO_NO_SCAN="1")
    p = run("doctor", py=venv_without_bio, env=env)
    assert p.returncode == 1 and "!! Biopython" in p.stdout
