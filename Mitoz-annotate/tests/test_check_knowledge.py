"""Tests for benchmark/check_knowledge.py (knowledge-base files must follow the template and cite sources).

Run: python3 -m pytest -q tests/test_check_knowledge.py
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "benchmark"))
import check_knowledge as ck  # noqa: E402

HEADINGS = ["## Version checked", "## What it does", "## Algorithm & key source locations", "## Defaults that matter",
            "## Known bugs / issues", "## Limitations by taxon", "## Relevance to our pipeline", "## References"]


def page(sections, tmp_path, name="tool.md"):
    text = "# Tool\n\n" + "\n\n".join(f"{h}\n{sections.get(h, '')}" for h in HEADINGS if h not in sections.get("_drop", []))
    p = tmp_path / name
    p.write_text(text + "\n", encoding="utf-8")
    return p


def ok_sections(**over):
    s = {"## Version checked": "- MitoZ 3.6, checked 2026-10-20 (https://github.com/linzhi2013/MitoZ)",
         "## What it does": "- annotates mitogenomes",
         "## Relevance to our pipeline": "- used for PCG and rRNA"}
    s.update(over)
    return s


def test_missing_heading(tmp_path):
    p = page(ok_sections(_drop=["## Known bugs / issues"]), tmp_path)
    assert "missing heading: ## Known bugs / issues" in ck.check_file(p)


def test_unsourced_bullet(tmp_path):
    p = page(ok_sections(**{"## Algorithm & key source locations": "- MiTFi is single-core"}), tmp_path)
    assert any("unsourced" in x for x in ck.check_file(p))


def test_sourced_bullets_ok(tmp_path):
    body = "\n".join(["- MiTFi runs with -cores 1 (mitoz/annotate/annotation.py:274)",
                      "- MITOS2 paper (doi:10.1093/nar/gkz833)",
                      "- crash on linear input (https://github.com/x/y/issues/12)",
                      "- may skip short tRNAs 未确证"])
    p = page(ok_sections(**{"## Algorithm & key source locations": body}), tmp_path)
    assert ck.check_file(p) == []


def test_version_line_required(tmp_path):
    p = page(ok_sections(**{"## Version checked": "- MitoZ 3.6 (https://github.com/linzhi2013/MitoZ)"}), tmp_path)
    assert any("checked YYYY-MM-DD" in x for x in ck.check_file(p))


def test_index_exempt(tmp_path):
    p = tmp_path / "INDEX.md"
    p.write_text("# Index\n\n- mitoz.md: MitoZ\n", encoding="utf-8")
    assert ck.check_file(p) == []
