#!/usr/bin/env python3
"""Run blast2metadata from a source checkout without installing (same as `python -m blast2metadata`)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from blast2metadata.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
