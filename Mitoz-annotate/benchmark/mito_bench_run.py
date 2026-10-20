"""Run annotation tools on the bare benchmark records and lift every result to NCBI coordinates.

Tools:
  ours    batch_mitoz.py run (current defaults) + mito_curate.py
  mitoz   MitoZ `annotate` unchanged (no post-processing)
  mitos2  MITOS2 runmitos, full mode (proteins, rRNA, tRNA), refseq89m
A record whose phylum MitoZ does not support (empty mitoz_clade) is 'unsupported' for ours and mitoz.

Per tool and record: <root>/runs/<tool>/<acc>/{final.gbf, status.json, log.txt}. final.gbf is in the
NCBI record's orientation and origin (lift_to_ncbi). Records with a status.json are skipped (resume).

Usage: python3 mito_bench_run.py --set sets/mito_bench_2026-10.tsv --root projects/mito_benchmark
         [--split dev] [--tools ours,mitoz,mitos2] [--jobs 4] [--threads 4] [--timeout 3600] [--acc A,B]
         [--recollect]
"""
import argparse
import csv
import json
import os
import re
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
PIPE = HERE.parent
sys.path.insert(0, str(PIPE))
import mito_env  # noqa: E402

TOOLS = ("ours", "mitoz", "mitos2")
ANCHOR = 50
MITOS_TRNA = re.compile(r"^trn([A-Z])(\d?)(?:_(\d+))?(?:-([a-z]))?\(([a-z-]{3})\)", re.I)
MITOS_GENE = re.compile(r"^(.+?)(?:_(\d+))?(?:-([a-z]))?$")       # nad5, nad5_1 (2nd copy), cox2-a (fragment)
MITOS_SKIP = re.compile(r"^(OH|OL)", re.I)


class LiftError(Exception):
    pass


# ── coordinates ──────────────────────────────────────────────────────────────
def _merge(parts: List[FeatureLocation]) -> List[FeatureLocation]:
    """Join consecutive parts that touch (an origin split undone by rotation)."""
    out: List[FeatureLocation] = []
    for p in parts:
        if out and out[-1].strand == p.strand:
            q = out[-1]
            if p.strand == -1 and int(p.end) == int(q.start):
                out[-1] = FeatureLocation(p.start, q.end, -1)
                continue
            if p.strand != -1 and int(q.end) == int(p.start):
                out[-1] = FeatureLocation(q.start, p.end, p.strand)
                continue
        out.append(p)
    return out


def _shift(part: FeatureLocation, off: int, n: int) -> List[FeatureLocation]:
    """Move one part by off on a circle of length n; split at the origin (biological order kept)."""
    s = (int(part.start) + off) % n
    e = s + int(part.end) - int(part.start)
    if e <= n:
        return [FeatureLocation(s, e, part.strand)]
    a, b = FeatureLocation(s, n, part.strand), FeatureLocation(0, e - n, part.strand)
    return [b, a] if part.strand == -1 else [a, b]


def lift_to_ncbi(pred: SeqRecord, ncbi_seq: Seq) -> SeqRecord:
    """Return pred re-expressed on ncbi_seq (same molecule, possibly rotated and/or reverse-complemented)."""
    n = len(ncbi_seq)
    pseq = Seq(str(pred.seq))
    if len(pseq) != n:
        raise LiftError(f"length {len(pseq)} != NCBI {n}")
    ring = str(ncbi_seq).upper() * 2
    for rc in (False, True):
        q = str(pseq.reverse_complement() if rc else pseq).upper()
        hits = [m.start() for m in re.finditer(f"(?={re.escape(q[:ANCHOR])})", ring[:n + ANCHOR - 1])]
        hits = [h for h in hits if ring[h:h + n] == q]
        if len(hits) == 1:
            off = hits[0]
            break
    else:
        raise LiftError("no unique match of the predicted sequence on the NCBI sequence")
    out = SeqRecord(ncbi_seq, id=pred.id, name=pred.name, description=pred.description)
    out.annotations = dict(pred.annotations)
    for f in pred.features:
        parts = list(f.location.parts)
        if rc:
            parts = [FeatureLocation(n - int(p.end), n - int(p.start), -(p.strand or 1)) for p in parts]
        new = _merge([x for p in parts for x in _shift(p, off, n)])
        loc = new[0] if len(new) == 1 else CompoundLocation(new)
        out.features.append(SeqFeature(loc, type=f.type, qualifiers=dict(f.qualifiers)))
    return out


# ── MITOS2 output ────────────────────────────────────────────────────────────
def mitos_to_record(result_dir: Path, seq: Seq, acc: str) -> SeqRecord:
    """GenBank-like record from MITOS2 result.bed (0-based). Copies (_0, _1) stay separate features of the
    same gene; fragments (-a, -b) are joined in letter order; start > end wraps the origin; OH/OL skipped."""
    n = len(seq)
    groups: Dict[tuple, List[tuple]] = {}
    for line in (result_dir / "result.bed").read_text().splitlines():
        c = line.split("\t")
        if len(c) < 6 or MITOS_SKIP.match(c[3]):
            continue
        strand = -1 if c[5].strip() == "-" else 1
        m = MITOS_TRNA.match(c[3])
        if m:
            key = ("tRNA", f"trn{m.group(1).upper()}{m.group(2)}", m.group(5).lower().strip("-"), m.group(3) or "")
            frag = m.group(4) or ""
        else:
            g = MITOS_GENE.match(c[3])
            base = g.group(1) if g else c[3]
            key = ("rRNA" if base.lower().startswith("rrn") else "CDS", base, "", (g.group(2) or "") if g else "")
            frag = (g.group(3) or "") if g else ""
        s0, e0 = int(c[1]), int(c[2])
        if s0 <= e0:
            parts = [FeatureLocation(s0, e0, strand)]
        else:                                    # circular feature across the origin
            a, b = FeatureLocation(s0, n, strand), FeatureLocation(0, e0, strand)
            parts = [b, a] if strand == -1 else [a, b]
        groups.setdefault(key, []).append((frag, parts))
    r = SeqRecord(seq, id=acc)
    r.annotations["molecule_type"] = "DNA"
    for (kind, name, ac, _), frags in groups.items():
        parts = [p for _, ps in sorted(frags, key=lambda x: x[0]) for p in ps]
        loc = parts[0] if len(parts) == 1 else CompoundLocation(parts)
        q = {"gene": [name]}
        if kind == "tRNA":
            q["anticodon"] = [ac]
        r.features.append(SeqFeature(loc, type="gene", qualifiers={"gene": [name]}))
        r.features.append(SeqFeature(loc, type=kind, qualifiers=q))
    return r


# ── running ──────────────────────────────────────────────────────────────────
def _env_path(exe: Optional[str]) -> Dict[str, str]:
    env = dict(os.environ)
    if exe:
        env["PATH"] = str(Path(exe).parent) + os.pathsep + env.get("PATH", "")
    return env


def build_cmd(tool: str, row: dict, bare: Path, out: Path, threads: int) -> Optional[List[List[str]]]:
    """Command sequence for one tool and record; None when the tool does not support the record's phylum."""
    acc, code = row["accession"], str(row["genetic_code"])
    if tool in ("ours", "mitoz") and not row["mitoz_clade"]:
        return None
    if tool == "ours":
        return [[sys.executable, str(PIPE / "batch_mitoz.py"), "run", "-i", str(bare), "-o", str(out / "raw"),
                 "--clade", row["mitoz_clade"], "--genetic_code", code, "--threads", str(threads), "--max_tasks", "1"],
                [sys.executable, str(PIPE / "mito_curate.py"), "-i", str(out / "raw" / "04.final_gb"),
                 "-o", str(out / "curated"), "--genetic_code", code, "--no_validate"]]
    fa = str(out / "in.fa")
    if tool == "mitoz":
        exe = mito_env.find("mitoz", None) or "mitoz"
        return [[exe, "annotate", "--workdir", str(out / "w"), "--outprefix", acc, "--thread_number", str(threads),
                 "--clade", row["mitoz_clade"], "--genetic_code", code, "--fastafiles", fa]]
    if tool == "mitos2":
        (out / "mitos").mkdir(parents=True, exist_ok=True)          # runmitos refuses a missing -o directory
        exe = mito_env.find("runmitos", None) or "runmitos"
        cmd = [exe, "-i", fa, "-c", code, "-o", str(out / "mitos"), "-r", "refseq89m",
               "--refdir", str(mito_env.mitos_refdir(None, "refseq89m") or ""), "--noplots"]
        return [cmd + (["--linear"] if row.get("topology") == "linear" else [])]
    raise ValueError(tool)


def _killpg(p: subprocess.Popen) -> None:
    """SIGTERM the process group, SIGKILL whatever is left after 5 s."""
    for sig, wait in ((signal.SIGTERM, 5), (signal.SIGKILL, 5)):
        try:
            os.killpg(p.pid, sig)
        except ProcessLookupError:
            return
        try:
            p.wait(wait)
            os.killpg(p.pid, 0)              # leader gone; children may still be in the group
        except (subprocess.TimeoutExpired, ProcessLookupError, PermissionError):
            pass
        time.sleep(0.2)


def run_steps(cmds: List[List[str]], out: Path, timeout: int) -> dict:
    """Run commands in order (each with its program's directory first on PATH); status dict."""
    out.mkdir(parents=True, exist_ok=True)
    st: dict = {"status": "ok", "seconds": 0.0, "step_seconds": {}, "error": ""}
    with (out / "log.txt").open("a") as log:
        for i, cmd in enumerate(cmds):
            t0 = time.time()
            log.write("$ " + " ".join(cmd) + "\n")
            log.flush()
            try:
                # own process group, so a timeout also stops the tools it started (MitoZ -> perl -> tblastn)
                p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.PIPE, text=True, errors="replace",
                                     start_new_session=True,
                                     env=_env_path(cmd[0] if os.sep in cmd[0] else None))
            except OSError as e:
                st.update(status="failed", error=f"step {i + 1}: {e}")
            else:
                try:
                    _, err = p.communicate(timeout=timeout)
                except subprocess.TimeoutExpired:
                    _killpg(p)
                    st.update(status="timeout", error=f"step {i + 1} exceeded {timeout} s")
                else:
                    log.write(err)
                    if p.returncode != 0:
                        st.update(status="failed", error=f"step {i + 1} exit {p.returncode}: {err.strip()[-500:]}")
            st["step_seconds"][Path(cmd[1] if cmd[0] == sys.executable else cmd[0]).name] = round(time.time() - t0, 1)
            st["seconds"] = round(st["seconds"] + time.time() - t0, 1)
            if st["status"] != "ok":
                break
    return st


def write_status(out: Path, st: dict) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / "status.json").write_text(json.dumps(st, indent=1))


def collect(tool: str, out: Path, ncbi: SeqRecord) -> SeqRecord:
    """The tool's annotation of this record, before lift-over."""
    if tool == "mitos2":
        return mitos_to_record(out / "mitos", Seq(str(ncbi.seq)), str(ncbi.id))
    pattern = {"ours": "curated/*.gbf", "mitoz": "w/**/*.gbf"}[tool]
    files = sorted(out.glob(pattern), key=lambda p: p.stat().st_size, reverse=True)
    if not files:
        raise LiftError(f"no GenBank output ({pattern})")
    return SeqIO.read(str(files[0]), "genbank")


def run_one(tool: str, row: dict, root: Path, threads: int, timeout: int, recollect: bool = False) -> dict:
    """Run one tool on one record (skipped when status.json exists). recollect: rebuild final.gbf from the
    kept tool output of an 'ok' run without running the tool again (after a parser fix)."""
    acc = row["accession"]
    out = root / "runs" / tool / acc
    old = json.loads((out / "status.json").read_text()) if (out / "status.json").exists() else None
    if old and not (recollect and old["status"] == "ok"):
        return old
    ncbi = SeqIO.read(str(root / "gb" / f"{acc}.gb"), "genbank")
    out.mkdir(parents=True, exist_ok=True)
    if old:
        st = {k: old[k] for k in ("status", "seconds", "step_seconds", "error")}
    else:
        (out / "in.fa").write_text(f">{acc}\n{ncbi.seq}\n")
        cmds = build_cmd(tool, row, root / "bare" / acc, out, threads)
        if cmds is None:
            st = {"status": "unsupported", "seconds": 0.0, "step_seconds": {}, "error": "phylum not supported by MitoZ"}
        else:
            st = run_steps(cmds, out, timeout)
    if st["status"] == "ok":
        try:
            lifted = lift_to_ncbi(collect(tool, out, ncbi), Seq(str(ncbi.seq)))
            SeqIO.write(lifted, str(out / "final.gbf"), "genbank")
        except Exception as e:                   # any collect / lift problem is this record's failure, not the run's
            st.update(status="failed", error=f"collect/lift: {type(e).__name__}: {e}")
    st = {"tool": tool, "accession": acc, **st}
    write_status(out, st)
    return st


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--set", required=True)
    p.add_argument("--root", required=True)
    p.add_argument("--split", default="dev", help="dev, test or all")
    p.add_argument("--tools", default=",".join(TOOLS))
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--timeout", type=int, default=3600)
    p.add_argument("--acc", help="only these accessions (comma-separated)")
    p.add_argument("--recollect", action="store_true",
                   help="rebuild final.gbf of finished runs from their kept tool output (no tool is re-run)")
    a = p.parse_args(argv)
    rows = list(csv.DictReader(open(a.set, encoding="utf-8"), delimiter="\t"))
    rows = [r for r in rows if a.split == "all" or r["split"] == a.split]
    if a.acc:
        rows = [r for r in rows if r["accession"] in a.acc.split(",")]
    tasks = [(t, r) for r in rows for t in a.tools.split(",")]
    root = Path(a.root)
    with ThreadPoolExecutor(a.jobs) as ex:
        for st in ex.map(lambda tr: run_one(tr[0], tr[1], root, a.threads, a.timeout, a.recollect), tasks):
            print(f"{st['tool']:7s} {st['accession']:14s} {st['status']:11s} {st['seconds']:8.1f}s {st['error'][:100]}",
                  flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
