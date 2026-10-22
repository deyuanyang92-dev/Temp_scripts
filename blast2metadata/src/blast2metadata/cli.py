"""Command-line interface."""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional


from . import __version__, ratelimit
from .doctor import run_doctor
from .engines import prepare_engine_db
from .fastautil import detect_sequence_type
from .fields import BLAST_REQUIRED_MODES
from .util import HAVE_PANDAS, HAVE_PSUTIL
from .blasttab import _sort_spec_arg
from .pipeline import add_header_to_table, find_input_files, preflight_checks, process_one_fasta, setup_output_directory
from .ratelimit import _make_blast_limiter, _make_entrez_limiter
from .smk import run_with_snakemake
from .util import ensure_dir, logger, resolve_total_threads

# ─────────────────────────────────────────────
# CLI (Fix-J: updated epilog examples)
# ─────────────────────────────────────────────
# Codex 修改(help-short): 新增短帮助文本生成函数，供 -h 和无参数错误提示复用。
# 理由: 日常运行只需要核心用法；完整高级参数仍由 argparse 的 --help 输出。
# 原代码摘要: 原来 -h/--help 都直接使用 argparse 自动生成的完整帮助。
def short_help_text(script_name: str) -> str:
    return f"""usage:
  {script_name} -i INPUT -o OUT [options]

常用一键运行:
  # 默认 remote / online BLAST
  {script_name} -i query.fasta -o OUT

  # 本地 BLAST 数据库
  {script_name} -i query.fasta -o OUT -bm local --blast_db_nt my_db

  # 本地库改用 MMseqs2 (核酸/蛋白) 或 DIAMOND (仅蛋白库, blastp/blastx); 库可直接给 FASTA
  {script_name} -i query.fasta -o OUT -bm local --engine mmseqs --blast_db refs.fasta
  {script_name} -i query.fasta -o OUT -bm local --engine diamond --blast_db nr.dmnd

  # 不使用 Snakemake，直接普通 Python 运行
  {script_name} -i query.fasta -o OUT --no_smk

  # 只从 accession list 获取 metadata
  {script_name} -i accessions.tsv -o OUT --mode metadata-only

  # 已有 BLAST TSV，只做 TopN
  {script_name} -i blast.tsv -o OUT --mode topN-only --topN 10

  # 给已有 BLAST 表加表头（Excel 查看）；列不是默认格式时用 --outfmt_blast 写明当时的 -outfmt
  {script_name} -i blast.tsv -o OUT --mode add-header --header_style both
  {script_name} -i blast.tsv -o OUT --mode add-header --outfmt_blast "6 std staxids"

核心参数:
  -i, --input FILE/DIR        输入 FASTA / BLAST TSV / accession list
  -o, --output DIR            输出目录，默认 BLAST_results
  --mode MODE                 默认 blast2topN2metadata
  -bm, --blast_method METHOD  online 或 local，默认 online
  --online_engine ENGINE      默认 blast_remote
  --blast_db_nt DB            本地核酸库
  --blast_db_nr DB            本地蛋白库
  --engine ENGINE             本地引擎: blast(默认)/diamond/mmseqs；--blast_db 可为 .dmnd、MMseqs2 库或 FASTA
  --doctor                    检查 blastn/diamond/mmseqs/snakemake 与 Python 依赖
  --topN N                    每个 query 保留前 N 个 hit，默认 10
  --sort_by SPEC              排序: bitscore(默认)/identity/coverage/evalue/pident_x_qcovs，
                              或自定义列顺序如 qcovs,pident
  --min_pident X / --min_qcovs X   TopN 前先过滤低相似度/低覆盖度命中
  --header_style STYLE        表头: code(pident) / label(英文说明) / both，默认 code
  --resume                    断点续跑
  --no_smk                    不用 Snakemake
  --email EMAIL               NCBI metadata/qblast 邮箱
  --no_ncbi_metadata          不联网获取 NCBI metadata
  --local_meta_tsv FILE       使用本地 metadata TSV

说明:
  默认就是 online/remote 流程；如果要本地库，加 -bm local。
  默认启用 Snakemake 调度；如果不想用，加 --no_smk。
  查看全部高级参数:
    {script_name} --help
"""


def parse_arguments(argv: Optional[List[str]] = None) -> argparse.Namespace:
    argv = list(sys.argv[1:] if argv is None else argv)
    script_name = "blast2metadata"

    # Codex 修改(help-short): -h 输出短帮助，--help 保留完整 argparse 帮助。
    # 理由: 用户日常只需要核心命令；高级参数列表很长，继续放在 --help。
    # 原代码摘要: 原来 argparse 默认让 -h 和 --help 显示同一份完整帮助。
    if "-h" in argv and "--help" not in argv:
        print(short_help_text(script_name))
        raise SystemExit(0)

    # Codex 修改(help-short): 无参数运行时显示短帮助后按缺少输入报错退出。
    # 理由: 比完整 argparse usage 更易读，同时保留 CLI 参数错误的退出码 2。
    # 原代码摘要: 原来无参数直接由 argparse 报 required -i/--input 错误。
    if not argv:
        print(short_help_text(script_name), file=sys.stderr)
        print("error: the following arguments are required: -i/--input", file=sys.stderr)
        raise SystemExit(2)

    p = argparse.ArgumentParser(
        prog="blast2metadata",
        description=(
            "Batch BLAST with TopN + metadata merge, local-meta enabled.\n\n"
            "Key notes for local NT-like DB:\n"
            "  1) local_nt_metadata.tsv usually uses accession as key\n"
            "  2) sequence_metadata.tsv usually uses pseudo_accession as key\n"
            "  3) default --local_meta_key auto detects the right column\n"
            "  4) Snakemake scheduling is ON by default for safe global resource control"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "Examples\n"
            "--------\n"
            "Pure local DB:\n"
            "  blast2metadata -i query.fasta -o OUT \\\n"
            "    --mode blast2topN2metadata -bm local --blast_db_nt my_local_db \\\n"
            "    --no_ncbi_metadata --local_meta_tsv ref_out/local_nt_metadata.tsv\n\n"
            "Mixed nt + local alias:\n"
            "  blast2metadata -i query.fasta -o OUT \\\n"
            "    --mode blast2topN2metadata -bm local --blast_db_nt nt_plus_local \\\n"
            "    --local_meta_tsv ref_out/local_nt_metadata.tsv --pseudo_acc_prefixes LOC\n\n"
            "Default global scheduling (auto CPU/RAM):\n"
            "  blast2metadata -i fasta_dir -o OUT\n\n"
            "Shared-server safe mode with busy-guard:\n"
            "  blast2metadata -i fasta_dir -o OUT \\\n"
            "    --threads 64 --smk_busy_guard --smk_auto_min_mem_mb_per_job 32000\n"
        ),
    )

    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("--doctor", action="store_true",
                   help="Report which external tools (BLAST+, DIAMOND, MMseqs2, snakemake) and Python "
                        "packages are available, then exit")
    p.add_argument("-i", "--input",  required=False, default="",
                   # Codex 修改: help 补充 topN-only 可直接接收 BLAST TSV。
                   # 理由: 修复后的行为覆盖旧脚本分步使用场景，帮助信息需要同步。
                   help="Input FASTA file/directory, BLAST TSV (topN-only), or accession list (metadata-only).")
    p.add_argument("-o", "--output", default="BLAST_results",
                   help="Output directory (default: BLAST_results)")
    p.add_argument("--resume", action="store_true", help="Resume from previous run")

    p.add_argument(
        "--mode",
        choices=["blast-only", "blast2topN", "blast2meta", "blast2topN2metadata",
                 "metadata-only", "topN-only", "add-header"],
        default="blast2topN2metadata",
        help="Operation mode. add-header: add a header row to existing BLAST tabular file(s) "
             "(outfmt 6/7/10/20) and write OUT/<name>_with_header.csv",
    )

    p.add_argument("--extension_fasta", default=".fasta,.fas,.faa,.fa,.fna",
                   help="Comma-separated FASTA extensions (default: .fasta,.fas,.faa,.fa,.fna)")
    p.add_argument("--suffix_fasta", default="", help="Optional suffix for FASTA files")

    # BLAST
    p.add_argument("-bm", "--blast_method", choices=["local", "online"], default="online",
                   help="BLAST method (default: online)")
    p.add_argument("--online_engine", choices=["blast_remote", "biopython"],
                   default="blast_remote",
                   help="Online engine: blast_remote (BLAST+ -remote) or biopython (qblast)")
    p.add_argument("--engine", choices=["blast", "diamond", "mmseqs"], default="blast",
                   help="Local search engine (-bm local): blast = BLAST+ (default); diamond = DIAMOND "
                        "(protein DB only: blastp/blastx, 100-10000x faster than BLAST for protein); "
                        "mmseqs = MMseqs2 (nucleotide or protein DB). DB from --blast_db: DIAMOND >= 2.1.14 "
                        "reads BLAST protein DBs directly (nr, swissprot; bare names are looked up in $BLASTDB) "
                        "or a .dmnd; MMseqs2 takes a DB from `mmseqs createdb` / `mmseqs databases NR|NT`; "
                        "both accept a FASTA, converted once into OUT/_engine_db/")
    p.add_argument("--query_gencode", type=int, default=None,
                   help="Genetic code for translated searches of nucleotide queries (blastx -query_gencode, "
                        "diamond --query-gencode, mmseqs --translation-table), e.g. 5 = invertebrate mito, "
                        "2 = vertebrate mito, 4 = mold/coelenterate mito. Default: tool default (1)")
    p.add_argument("--diamond_exe", default="diamond", help="DIAMOND executable (default: diamond)")
    p.add_argument("--mmseqs_exe",  default="mmseqs",  help="MMseqs2 executable (default: mmseqs)")
    p.add_argument("--blast_exe",    default="auto",
                   help="BLAST executable/program name (auto → blastn/blastp)")
    p.add_argument("--blast_db",     default="auto",
                   help="BLAST database name/path (auto → nt/nr)")
    p.add_argument("--blast_db_nt",  default="nt",
                   help="Local nucleotide BLAST database path/name (default: nt)")
    p.add_argument("--blast_db_nr",  default="nr",
                   help="Local protein BLAST database path/name (default: nr)")
    p.add_argument("--evalue",           default=1e-5,  type=float, help="E-value threshold")
    p.add_argument("--max_target_seqs",  default=100,   type=_positive_int,
                   help="Max target seqs per query (default: 100)")
    p.add_argument("--outfmt_blast",     default="6",
                   help='Output format: "6" or "7" uses built-in 14 columns; '
                        'custom: "6 qseqid sseqid …" must include all 14 standard fields '
                        '(keyword "std" is expanded). In add-header mode: the -outfmt you ran BLAST with, '
                        'e.g. "6 std staxids"')
    p.add_argument("--header_style", choices=["code", "label", "both"], default="code",
                   help='Header text of *_results.csv / *_with_header.csv: code = "pident"; '
                        'label = "Percentage of identical matches"; both = "pident (Percentage of '
                        'identical matches)" (default: code)')
    p.add_argument("--no_blast_header", action="store_true",
                   help="Do not write a header row in blast_results/<name>_blast.tsv")
    p.add_argument("-epb", "--extra_para2_blast", default="",
                   help="Extra parameters appended to BLAST command (do NOT include -outfmt)")

    # Resources
    p.add_argument("--threads", default=0, type=int,
                   help="Total CPU thread ceiling. 0 = auto-detect host CPUs")
    p.add_argument("--parallel_jobs", default=None, type=int,
                   help="Concurrent BLAST jobs inside one script run")
    p.add_argument("--local_threads_per_job", default=None, type=int,
                   help="Override per-job -num_threads for local BLAST")
    p.add_argument("--batch_size",            default=None, type=_positive_int,
                   help="Sequences per search batch (default: 5000 for local BLAST+/DIAMOND/MMseqs2, 50 online). "
                        "Each local batch is one full database pass, so fewer, larger batches are faster")
    p.add_argument("--online_max_concurrent", default=1,  type=int,
                   help="Max concurrent online BLAST batches (default: 1)")
    # 修改 (2026-10-23): 默认 12 → 0。提交间隔已由 BLAST 限速器保证 ≥10 s（NCBI 规定），
    # 这段固定等待是重复的；原来 qblast 逐条提交时每条都要多等 12 s。
    p.add_argument("--online_delay_sec",      default=0, type=int,
                   help="Extra delay in seconds after each online batch (default: 0; submissions are "
                        "already spaced >=10 s apart as NCBI requires)")
    p.add_argument("--retries",    default=3,  type=_positive_int,
                   help="Retries for failed BLAST/metadata fetch (default: 3)")
    p.add_argument("--blast_timeout", default=0, type=_nonnegative_float,
                   help="Timeout in seconds per external search attempt; 0 disables it (not qblast)")
    p.add_argument("--retry_sleep", default=10, type=int,
                   help="Initial retry sleep in seconds (default: 10)")

    # TopN
    p.add_argument("--topN",    default=10, type=_positive_int,
                   help="Top N hits per query (default: 10)")
    # 修改 (2026-10-10): 默认 identity -> bitscore。理由: identity 优先会让 30 bp/100% 的短命中
    # 排在全长 99% 命中之前；bitscore 同时反映长度和相似度。
    # 修改 (2026-10-11): 可用预设或自定义多列排序，并加阈值过滤。
    p.add_argument("--sort_by", type=_sort_spec_arg, default="bitscore",
                   help="TopN ranking. Preset: bitscore (default), identity, coverage, evalue, "
                        "pident_x_qcovs (= pident*qcovs/100, ranking heuristic); "
                        "or comma-separated columns ranked in order, e.g. 'qcovs,pident' "
                        "(evalue/mismatch/gapopen/gaps: smaller is better; others: larger is better)")
    p.add_argument("--min_pident", type=float, default=0.0,
                   help="Drop hits with pident below this before TopN (default: 0 = off)")
    p.add_argument("--min_qcovs", type=float, default=0.0,
                   help="Drop hits with qcovs below this before TopN (default: 0 = off)")

    # Entrez
    p.add_argument("--email",       default="you@example.com",
                   help="Email address for NCBI Entrez / qblast")
    p.add_argument("--api_key",     default="",
                   help="NCBI API key (enables 10 req/sec for Entrez E-utilities)")
    p.add_argument("--entrez_tool", default="blast2metadata",
                   help="Entrez.tool string")
    p.add_argument("--entrez_delay_sec", default=1.0, type=float,
                   help="Extra delay between Entrez batches in seconds (default: 1.0)")

    p.add_argument("-mty", "--metadata_type", choices=["all", "short"], default="all",
                   help="Metadata verbosity: all or short (default: all)")
    p.add_argument("-sm",  "--shorter_metadata",
                   default="organism,reference_titles,reference_authors,reference_journals",
                   help="Fields included when --metadata_type short")

    # Local meta
    p.add_argument("--no_ncbi_metadata", action="store_true",
                   help="Skip Entrez efetch entirely")
    p.add_argument("--local_meta_tsv", default="",
                   help="Path to local metadata TSV (e.g. local_nt_metadata.tsv)")
    p.add_argument("--local_meta_key", default="auto",
                   help="Key column in local meta TSV. auto → accession > pseudo_accession > first col")
    p.add_argument("--local_meta_mode", choices=["map_to_ncbi", "prefix_only"],
                   default="map_to_ncbi",
                   help="How to merge local meta: map_to_ncbi or prefix_only")
    p.add_argument("--local_meta_prefix", default="local_",
                   help="Column prefix when --local_meta_mode=prefix_only")
    p.add_argument("--pseudo_acc_prefixes", default="LOC",
                   help="Comma-separated accession prefixes treated as local/pseudo")

    # Snakemake, default ON
    p.add_argument("--smk", "--snakemake", dest="use_snakemake", action="store_true",
                   help="Enable built-in Snakemake scheduling (default: ON)")
    p.add_argument("--no_smk", dest="use_snakemake", action="store_false",
                   help="Disable Snakemake and run plain Python loop")

    p.add_argument("--smk_cores", type=int, default=0,
                   help="Global Snakemake cores ceiling. 0 = auto")
    p.add_argument("--smk_jobs",  type=int, default=0,
                   help="Global concurrent Snakemake jobs. 0 = auto")
    p.add_argument("--smk_threads_per_job", type=int, default=0,
                   help="Threads per rule. 0 = auto")
    p.add_argument("--smk_total_mem_mb",  type=int, default=0,
                   help="Global mem_mb ceiling. 0 = auto")
    p.add_argument("--smk_mem_mb_per_job", type=int, default=0,
                   help="Per-job mem_mb. 0 = auto")

    p.add_argument("--smk_auto_mem", dest="smk_auto_mem", action="store_true",
                   help="Auto-detect RAM for Snakemake mem_mb budget (default: ON)")
    p.add_argument("--no_smk_auto_mem", dest="smk_auto_mem", action="store_false",
                   help="Disable automatic RAM detection")
    p.add_argument("--smk_mem_fraction", type=float, default=0.90,
                   help="Fraction of total RAM considered (default: 0.90)")
    p.add_argument("--smk_max_available_mem_fraction", type=float, default=0.85,
                   help="Max fraction of currently available RAM exposed (default: 0.85)")
    p.add_argument("--smk_auto_min_threads_per_job", type=int, default=4,
                   help="Min threads kept per job in auto mode (default: 4)")
    p.add_argument("--smk_auto_min_mem_mb_per_job", type=int, default=24576,
                   help="Min memory per job in auto mode in MB (default: 24576)")
    p.add_argument("--smk_auto_max_jobs", type=int, default=8,
                   help="Auto mode upper bound for concurrent jobs (default: 8)")

    p.add_argument("--smk_total_job_slots",  type=int, default=0,
                   help="Abstract job_slots resource total. 0 = smk_jobs")
    p.add_argument("--smk_job_slots_per_job", type=int, default=1,
                   help="job_slots consumed per job (default: 1)")
    p.add_argument("--smk_latency_wait",  type=int, default=60,
                   help="Snakemake --latency-wait in seconds (default: 60)")
    p.add_argument("--smk_restart_times", type=int, default=1,
                   help="Snakemake restart/retry count for failed rules (default: 1)")

    p.add_argument("--smk_busy_guard", dest="smk_busy_guard", action="store_true",
                   help="Auto-throttle CPU/memory when server is busy (default: ON)")
    p.add_argument("--no_smk_busy_guard", dest="smk_busy_guard", action="store_false",
                   help="Disable busy-server throttling")
    p.add_argument("--smk_busy_guard_min_cores",  type=int,   default=1,
                   help="Busy guard: minimum global cores to keep (default: 1)")
    p.add_argument("--smk_busy_guard_min_mem_mb", type=int,   default=4096,
                   help="Busy guard: minimum global memory in MB to keep (default: 4096)")
    p.add_argument("--smk_idle_cpu_fraction",       type=float, default=0.90,
                   help="Idle server: max fraction of CPU to expose (default: 0.90)")
    p.add_argument("--smk_idle_avail_mem_fraction", type=float, default=0.85,
                   help="Idle server: max fraction of available RAM to expose (default: 0.85)")

    p.add_argument("--smk_keep_going", dest="smk_keep_going", action="store_true",
                   help="Snakemake --keep-going (default: ON)")
    p.add_argument("--no_smk_keep_going", dest="smk_keep_going", action="store_false",
                   help="Disable --keep-going")
    p.add_argument("--smk_rerun_incomplete", dest="smk_rerun_incomplete", action="store_true",
                   help="Snakemake --rerun-incomplete (default: ON)")
    p.add_argument("--no_smk_rerun_incomplete", dest="smk_rerun_incomplete", action="store_false",
                   help="Disable --rerun-incomplete")
    p.add_argument("--smk_printshellcmds", dest="smk_printshellcmds", action="store_true",
                   help="Snakemake --printshellcmds (default: ON)")
    p.add_argument("--no_smk_printshellcmds", dest="smk_printshellcmds", action="store_false",
                   help="Disable --printshellcmds")
    p.add_argument("--smk_unlock_on_lock_exception",
                   dest="smk_unlock_on_lock_exception", action="store_true",
                   help="Auto-run snakemake --unlock on stale lock (default: ON)")
    p.add_argument("--no_smk_unlock_on_lock_exception",
                   dest="smk_unlock_on_lock_exception", action="store_false",
                   help="Disable automatic --unlock recovery")
    p.add_argument("--smk_dry_run",  action="store_true", help="Snakemake --dry-run")
    p.add_argument("--smk_profile",  default=None,  help="Snakemake --profile path")
    p.add_argument("--smk_extra",    default="",    help="Extra snakemake CLI arguments")

    p.set_defaults(
        use_snakemake=True,
        smk_auto_mem=True,
        smk_keep_going=True,
        smk_rerun_incomplete=True,
        smk_printshellcmds=True,
        smk_unlock_on_lock_exception=True,
        smk_busy_guard=True,
    )

    # Hidden: set by per-stem and Snakemake child invocations
    p.add_argument("--out_stem",  default="", help=argparse.SUPPRESS)

    # Fix-1: child processes load the full serialised args from this JSON
    p.add_argument("--args_json", default="", help=argparse.SUPPRESS)

    # Fix-A/G: injected by parent into child_args.json; child uses it to scale
    # the per-process NCBI rate-limit share.  Never set by user directly.
    p.add_argument("--_smk_active_jobs", default=1, type=int, help=argparse.SUPPRESS)
    p.add_argument("--_engine_db", default="", help=argparse.SUPPRESS)        # set by prepare_engine_db
    p.add_argument("--_engine_db_type", default="", help=argparse.SUPPRESS)

    args = p.parse_args(argv)
    if not args.doctor and not args.input:
        p.error("the following arguments are required: -i/--input")
    return args


def _nonnegative_float(s: str) -> float:
    import math
    v = float(s)
    if not math.isfinite(v) or v < 0:
        raise argparse.ArgumentTypeError("must be finite and >= 0")
    return v


def _positive_int(s: str) -> int:
    try:
        v = int(s)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected an integer, got {s!r}")
    if v <= 0:
        raise argparse.ArgumentTypeError(f"must be >= 1, got {v}")
    return v


# ─────────────────────────────────────────────
# main (Fix-G: split ENTREZ_LIMITER and BLAST_LIMITER)
# ─────────────────────────────────────────────
def _absolutize_local_paths(args: argparse.Namespace) -> None:
    """Make local DB/metadata/exe paths absolute before any child process changes cwd.

    修改 (2026-10-10): Snakemake 子进程 cwd 为 OUT/_snakemake_workflow，原代码传相对路径
    （如 --blast_db_nt ../db/mito）时子任务报 "No alias or index file found"。
    远程库名（nt/nr、BLASTDB 中的库名）不含路径分隔符且当前目录无同名文件，保持原样。
    """
    for attr in ("blast_db", "blast_db_nt", "blast_db_nr"):
        v = str(getattr(args, attr, "") or "")
        if not v or v == "auto":
            continue
        local = os.sep in v or Path(v).exists() or any(
            Path(v + ext).exists() for ext in (".nal", ".pal", ".ndb", ".pdb", ".nin", ".pin", ".dmnd", ".dbtype"))
        if local:
            setattr(args, attr, str(Path(v).expanduser().resolve()))
    if getattr(args, "local_meta_tsv", ""):
        args.local_meta_tsv = str(Path(args.local_meta_tsv).expanduser().resolve())
    for attr in ("blast_exe", "diamond_exe", "mmseqs_exe"):
        if os.sep in str(getattr(args, attr, "")):
            setattr(args, attr, str(Path(getattr(args, attr)).expanduser().resolve()))


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_arguments(argv)
    if args.doctor:
        return run_doctor({"diamond": args.diamond_exe, "mmseqs": args.mmseqs_exe})

    # Fix-1: if child process, load base args from JSON then apply CLI overrides
    if args.args_json:
        json_path = Path(args.args_json)
        if not json_path.exists():
            logger.error(f"--args_json file not found: {json_path}")
            return 1
        with open(json_path, "r", encoding="utf-8") as fh:
            stored: Dict[str, object] = json.load(fh)

        # Fields the child always overrides from its own CLI
        child_override_keys = {"input", "output", "out_stem", "resume", "args_json", "use_snakemake"}
        child_overrides = {k: getattr(args, k) for k in child_override_keys if hasattr(args, k)}

        # Start from stored, apply child-specific overrides
        for k, v in stored.items():
            if k not in child_override_keys:
                setattr(args, k, v)
        for k, v in child_overrides.items():
            setattr(args, k, v)

        # Child jobs never recurse into Snakemake
        args.use_snakemake = False

    _absolutize_local_paths(args)
    for name in ("min_pident", "min_qcovs"):
        if not 0 <= getattr(args, name) <= 100:
            raise SystemExit(f"ERROR: --{name} must be between 0 and 100")
    if args.topN > args.max_target_seqs and args.mode in {"blast2topN", "blast2topN2metadata"}:
        logger.warning("[TOPN] topN exceeds max_target_seqs; the search may return too few candidates")
    args.threads = resolve_total_threads(args)
    if args.batch_size is None:
        # 修改 (2026-10-23): 本地 BLAST+ 默认 500 → 5000（与 DIAMOND/MMseqs2 一致）。
        # 理由: BLAST 每个批次要把整个数据库扫一遍，耗时主要在扫库而不是 query 数。
        #       实测 (blastn 2.12, nt 前 10 卷 28 GB 已在缓存, 8 线程, 600 bp COI):
        #       6 条 13.7 s, 96 条 20.0 s, 960 条 66.9 s —— 每批约 13 s 固定扫库开销;
        #       960 条切成 20 批 × 4 进程 211.9 s，1 批 66.9 s，HSP 结果完全相同。
        #       完整 nt 的 .nsq 有 800 GB，超过本机 503 GB 内存，每多一批就多读一遍磁盘。
        #       megablast 内部每次最多拼接约 5 Mbp query，5000 × 600 bp 仍在一遍之内；
        #       更长的 query 由 BLAST 自己内部分遍，不会因批大而变慢。
        #       代价: 单批失败/中断时 --resume 要重跑更多序列，可用 --batch_size 调小。
        args.batch_size = 5000 if args.blast_method == "local" else 50
    out_dir = ensure_dir(Path(args.output).resolve())
    args.output = str(out_dir)

    # Fix-G: reinitialise BOTH limiters now that we know the real args.
    # ENTREZ_LIMITER:
    #   • api_key  → NCBI allows 10 req/sec (instead of 3) for E-utilities
    #   • _smk_active_jobs → injected by the parent process into child_args.json;
    #     each child divides the global budget evenly so 10 parallel jobs never
    #     collectively exceed 3 (or 10) req/sec.
    # BLAST_LIMITER:
    #   • NCBI BLAST URL API: at most 1 request / 10 s, no API key support
    #   • _smk_active_jobs → divides this limit across parallel jobs
    active_jobs    = int(getattr(args, "_smk_active_jobs", 1))
    
    ratelimit.ENTREZ_LIMITER = _make_entrez_limiter(api_key=getattr(args, "api_key", ""),
                                          active_jobs=active_jobs)
    ratelimit.BLAST_LIMITER = _make_blast_limiter(active_jobs=active_jobs)
    
    _entrez_max_reqs = 10 if getattr(args, "api_key", "") else 3

    fh = logging.FileHandler(out_dir / "batch_blast.log", encoding="utf-8")
    fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(fh)

    exts  = [x.strip() for x in args.extension_fasta.split(",") if x.strip()]
    files = find_input_files(args.input, exts, args.suffix_fasta, args.mode)
    if not files:
        logger.error("No input files found.")
        return 1

    if args.mode == "add-header":
        bad = sum(not add_header_to_table(Path(fp), args, out_dir) for fp in files)
        return 1 if bad else 0

    preflight_checks(args, files, out_dir)

    logger.info(f"[FILES]          {len(files)}")
    logger.info(f"[PANDAS]         {HAVE_PANDAS}")
    logger.info(f"[PSUTIL]         {HAVE_PSUTIL}")
    logger.info(f"[THREAD_CEILING] {args.threads}")
    logger.info(f"[SCHEDULER]      {'snakemake' if args.use_snakemake else 'python'}")
    logger.info(
        f"[ENTREZ_RATE_LIMIT] {_entrez_max_reqs} calls / "
        f"{1.05 * active_jobs:.2f}s period "
        f"(api_key={'yes' if getattr(args,'api_key','') else 'no'}, "
        f"active_jobs={active_jobs})"
    )
    logger.info(
        f"[BLAST_RATE_LIMIT] 1 call / "
        f"{10.1 * active_jobs:.2f}s period "
        f"(strict NCBI BLAST limit, no API key support, "
        f"active_jobs={active_jobs})"
    )

    if args.use_snakemake:
        # 修改 (2026-10-22): DIAMOND/MMseqs2 数据库在父进程先建一次，子任务直接复用。
        if args.mode in BLAST_REQUIRED_MODES and args.blast_method == "local" and args.engine != "blast":
            prepare_engine_db(args, detect_sequence_type(files[0]), out_dir)
        return int(run_with_snakemake(args, files, out_dir))

    dirs = setup_output_directory(out_dir)
    # Codex 修改: 记录普通 Python 模式下的文件级失败数。
    # 理由: 原代码只记录异常并继续，最后仍返回 0，自动化流程会把失败任务误判为成功。
    # 原代码:
    #   for fp in files:
    #       try:
    #           process_one_fasta(fp, args, dirs)
    #       except ...
    failed = 0
    for fp in files:
        try:
            if not process_one_fasta(fp, args, dirs):
                failed += 1
        except KeyboardInterrupt:
            logger.error("Interrupted by user.")
            return 130
        except Exception as e:
            failed += 1
            logger.error(f"[FILE_FAIL] {fp}: {e}", exc_info=True)

    if failed:
        logger.error(f"Processing completed with {failed} failed file(s).")
        return 1

    logger.info("All processing completed.")
    return 0
