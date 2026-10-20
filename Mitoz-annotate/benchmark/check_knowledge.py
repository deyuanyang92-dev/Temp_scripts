"""Check knowledge-base Markdown files: template headings present, every claim bullet cites a source.

A bullet is sourced when it contains (path/file.py:123), (doi:...), (https://...), or is marked 未确证.
Bullets under "What it does" and "Relevance to our pipeline" need no source. INDEX.md and
benchmark_findings.md are free-form. Template headings apply to tools/*.md only ("Version checked" must say
"checked YYYY-MM-DD"); topics/*.md may use any headings but every bullet still needs a source.

Usage: python3 check_knowledge.py <knowledge_dir_or_file> ...   (exit 1 on any problem)
"""
import re
import sys
from pathlib import Path
from typing import Dict, List

HEADINGS = ["## Version checked", "## What it does", "## Algorithm & key source locations", "## Defaults that matter",
            "## Known bugs / issues", "## Limitations by taxon", "## Relevance to our pipeline", "## References"]
NO_SOURCE = {"## What it does", "## Relevance to our pipeline"}
FREE_FORM = {"INDEX.md", "benchmark_findings.md", "_TEMPLATE.md"}
SOURCE = re.compile(r"\(([\w./+-]+:\d+(-\d+)?|doi:\S+|https?://\S+)\)|未确证")
CHECKED = re.compile(r"checked \d{4}-\d{2}-\d{2}")


def check_file(path: Path) -> List[str]:
    """Problems found in one file (empty list = OK)."""
    if path.name in FREE_FORM:
        return []
    out = []
    section = ""
    sections: Dict[str, List[str]] = {}
    for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if line.startswith("## "):
            section = line.strip()
            sections[section] = []
            continue
        if section:
            sections[section].append(line)
        if line.lstrip().startswith("- ") and section and section not in NO_SOURCE and not SOURCE.search(line):
            out.append(f"{path.name}:{n}: unsourced bullet in {section}: {line.strip()[:80]}")
    if path.parent.name != "tools":                  # topics/: free headings, sources still required
        return out
    for h in HEADINGS:
        if h not in sections:
            out.append(f"missing heading: {h}")
    if "## Version checked" in sections and not CHECKED.search("\n".join(sections["## Version checked"])):
        out.append(f"{path.name}: ## Version checked needs 'checked YYYY-MM-DD'")
    return out


def main(argv: List[str]) -> int:
    files = []
    for a in argv:
        p = Path(a)
        files += sorted(p.rglob("*.md")) if p.is_dir() else [p]
    problems = [x for f in files for x in check_file(f)]
    for x in problems:
        print(x)
    print(f"{len(files)} files, {len(problems)} problems")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
