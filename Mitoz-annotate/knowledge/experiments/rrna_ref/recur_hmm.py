"""Add mito_rrna_hmm (LIB env) then re-curate kept raw outputs of a tool with the current mito_curate.py into <root2>/runs/<newtool>/ (offline A/B)."""
import csv, json, subprocess, sys, os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
PIPE = Path("/mnt/n/yzz-分类工作流/temp_scripts/Mitoz-annotate")
sys.path.insert(0, str(PIPE / "benchmark"))
import mito_bench_run as mr
sys.path.insert(0, str(PIPE))
import mito_rrna_hmm, shutil
LIB = Path(os.environ["LIB"])
REP = []
R = Path("/mnt/n/yzz-分类工作流/projects/mito_benchmark")
src_tool, root2, new = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
EXTRA = sys.argv[4:]
rows = {r["accession"]: r for r in csv.DictReader(open(PIPE / "benchmark/sets/mito_bench_2026-10.tsv"), delimiter="\t")}
for d in ("gb", "bare"):
    if not (root2 / d).exists(): (root2 / d).symlink_to(R / d)
(root2 / "runs").mkdir(parents=True, exist_ok=True)
for t in ("ours", "mitoz", "mitos2", "mitom"):
    if t != new and not (root2 / "runs" / t).exists() and (R / "runs" / t).exists(): (root2 / "runs" / t).symlink_to(R / "runs" / t)
def one(st):
    a = st.parent.name; o = json.load(open(st))
    out = root2 / "runs" / new / a; out.mkdir(parents=True, exist_ok=True)
    if o["status"] == "ok":
        fin = out / "04.final_gb"
        if fin.exists(): shutil.rmtree(fin)
        shutil.copytree(st.parent / "raw" / "04.final_gb", fin)
        REP.extend(dict(r, acc=a) for r in mito_rrna_hmm.search_dir(fin, LIB))
        p = subprocess.run([sys.executable, "-W", "ignore", str(PIPE / "mito_curate.py"), "-i", str(fin),
                            "-o", str(out / "curated"), "--genetic_code", rows[a]["genetic_code"], "--no_validate", *EXTRA], capture_output=True, text=True)
        if p.returncode: o = dict(o, status="failed", error=p.stderr[-200:])
        else:
            ncbi = SeqIO.read(str(R / "gb" / f"{a}.gb"), "genbank")
            SeqIO.write(mr.lift_to_ncbi(mr.collect("ours", out, ncbi), Seq(str(ncbi.seq))), str(out / "final.gbf"), "genbank")
    mr.write_status(out, dict(o, tool=new))
with ThreadPoolExecutor(6) as ex: list(ex.map(one, sorted((R / "runs" / src_tool).glob("*/status.json"))))
with open(root2 / f"{new}_rrna_hmm_report.tsv", "w") as h:
    h.write("acc\trecord\tphylum\tgene\tlocation\tstatus\n")
    for r in REP: h.write("\t".join(r[k] for k in ("acc", "record", "phylum", "gene", "location", "status")) + "\n")
print("done")
