"""The `mitom` command: one entry point for every MitoM step.

Run: python3 -m pytest -q tests/test_mitom_cli.py
"""
import subprocess
import sys
from pathlib import Path

MITOM = Path(__file__).resolve().parent.parent / "mitom"


def run(*a):
    return subprocess.run([sys.executable, str(MITOM), *a], capture_output=True, text=True, timeout=120)


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
