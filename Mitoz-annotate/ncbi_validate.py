#!/usr/bin/env python3
"""ncbi_validate.py - check GenBank annotation against NCBI's own submission validator.

Each record is converted to the NCBI submission inputs (.fsa with source modifiers and
5-column feature table .tbl) and validated with table2asn -V vb, the program NCBI uses
when it processes a submission. Messages are summarised per severity and code.
Submitter/publication messages (no template) are ignored: only sequence and feature
problems (SEQ_FEAT, SEQ_INST, SEQ_DESCR on the record) count.

Usage:
  python3 ncbi_validate.py -i gb_file_or_dir -o out_dir [--table2asn path] [--focus intron]
Outputs: out_dir/<record>/{.fsa,.tbl,.val}, out_dir/validation.tsv (one row per message),
and a summary on stdout.
"""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO
from Bio.SeqFeature import AfterPosition, BeforePosition

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402  (finds external programs: config, PATH, conda environments)

GB_EXT = {".gb", ".gbk", ".gbf", ".genbank", ".gbff"}
SKIP_QUALS = {"translation", "db_xref", "protein_id", "locus_tag_old"}
IGNORE_CODES = re.compile(r"SEQ_DESCR\.(NoPubFound|NoSourceDescriptor|MissingLineage|NoTaxonID|BadOrgMod)|"
                          r"GENERIC\.MissingPubRequirement|SEQ_PKG\.|SEQ_DESCR\.TitleNotAppropriateForSet")
MSG_RE = re.compile(r"^(INFO|WARNING|ERROR|REJECT|FATAL):\s+valid\s+\[([A-Z_]+\.[A-Za-z0-9_]+)\]\s+(.*)$", re.I)


def feature_lines(rec) -> List[str]:
    out = [f">Feature {rec.id}"]
    for f in rec.features:
        if f.type == "source":
            continue
        parts = list(f.location.parts)               # Biopython keeps minus-strand parts in transcript order
        for i, p in enumerate(parts):
            s, e = int(p.start) + 1, int(p.end)
            lo = "<" if isinstance(p.start, BeforePosition) else ""
            hi = ">" if isinstance(p.end, AfterPosition) else ""
            if p.strand == -1:
                a, b = f"{hi}{e}", f"{lo}{s}"
            else:
                a, b = f"{lo}{s}", f"{hi}{e}"
            out.append(f"{a}\t{b}" + (f"\t{f.type}" if i == 0 else ""))
        for k, vals in f.qualifiers.items():
            if k in SKIP_QUALS:
                continue
            for v in vals:
                out.append(f"\t\t\t{k}\t{' '.join(str(v).split())}")
    return out


def fasta_defline(rec, gencode: Optional[int]) -> str:
    org = rec.annotations.get("organism", "unknown organism")
    topo = rec.annotations.get("topology", "linear")
    mods = [f"[organism={org}]", "[location=mitochondrion]", f"[topology={topo}]"]
    if topo == "circular":                       # else SEQ_INST.CompleteCircleProblem (table2asn 2026-10-21)
        mods.append("[completeness=complete]")
    if gencode:
        mods.append(f"[mgcode={gencode}]")
    return f">{rec.id} " + " ".join(mods)


def record_gencode(rec) -> Optional[int]:
    codes = [int(f.qualifiers["transl_table"][0]) for f in rec.features
             if f.type == "CDS" and f.qualifiers.get("transl_table", [""])[0].isdigit()]
    return max(set(codes), key=codes.count) if codes else None


def our_intervals(rec, tag: str) -> List[str]:
    """'a-b' strings of every interval of features whose /note contains tag."""
    out = []
    for f in rec.features:
        if any(tag in n for n in f.qualifiers.get("note", [])):
            out += [f"{int(p.start) + 1}-{int(p.end)}" for p in f.location.parts]
    return out


def validate_record(rec, out: Path, table2asn: str, sbt: Optional[str] = None) -> List[Dict[str, str]]:
    """Write rec as NCBI submission inputs (rec.fsa + rec.tbl) and run table2asn on them (rec.sqn, rec.val);
    sbt = submitter template from https://submit.ncbi.nlm.nih.gov/genbank/template/submission/ (optional)."""
    d = out / re.sub(r"[^\w.]+", "_", rec.id)
    d.mkdir(parents=True, exist_ok=True)
    base = d / "rec"
    (d / "rec.fsa").write_text(fasta_defline(rec, record_gencode(rec)) + "\n" + str(rec.seq) + "\n")
    (d / "rec.tbl").write_text("\n".join(feature_lines(rec)) + "\n")
    subprocess.run([table2asn, "-i", str(d / "rec.fsa"), "-V", "vb", "-o", str(base) + ".sqn",
                    *(["-t", sbt] if sbt else [])],
                   capture_output=True, text=True, check=False)
    rows = []
    val = d / "rec.val"
    for line in val.read_text(errors="ignore").splitlines() if val.exists() else []:
        m = MSG_RE.match(line.strip())
        if not m or IGNORE_CODES.search(m.group(2)):
            continue
        mine = our_intervals(rec, "mito_intron")
        msg = m.group(3)
        rows.append({"record": rec.id, "severity": m.group(1).upper(), "code": m.group(2),
                     "ours": "yes" if any(re.search(rf"\b{re.escape(iv)}\b", msg) for iv in mine) else "no",
                     "message": msg[:300]})
    return rows


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True)
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--table2asn", default=None, help="table2asn executable (default: found by mito_env.py)")
    p.add_argument("--focus", default="", help="Only count messages whose text matches this regex (e.g. intron|CDS)")
    p.add_argument("--only_ours", action="store_true",
                   help="Only count messages on features written by mito_intron.py")
    a = p.parse_args(argv)
    exe = mito_env.find("table2asn", a.table2asn)
    if not exe:
        print("[ncbi_validate] table2asn not found (conda install -c bioconda table2asn)", file=sys.stderr)
        return 2
    src = Path(a.input)
    files = sorted(x for x in src.rglob("*") if x.suffix.lower() in GB_EXT) if src.is_dir() else [src]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, str]] = []
    n = 0
    for f in files:
        for rec in SeqIO.parse(str(f), "genbank"):
            n += 1
            rows += validate_record(rec, out, exe)
    with open(out / "validation.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["record", "severity", "code", "ours", "message"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    focus = re.compile(a.focus, re.I) if a.focus else None
    sel = [r for r in rows if (not focus or focus.search(r["message"])) and (not a.only_ours or r["ours"] == "yes")]
    sev = Counter(r["severity"] for r in sel)
    print(f"[ncbi_validate] {n} record(s); messages: " + ", ".join(f"{k} {v}" for k, v in sev.most_common()))
    for (s, c), k in Counter((r["severity"], r["code"]) for r in sel if r["severity"] in ("ERROR", "REJECT", "FATAL")).most_common(15):
        print(f"   {s:6s} {c:40s} {k}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
