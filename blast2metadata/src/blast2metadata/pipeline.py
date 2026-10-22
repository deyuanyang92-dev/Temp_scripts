"""Per-file pipeline: input discovery, preflight, BLAST -> TopN -> metadata, add-header."""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import glob
import json
import os
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple


from .util import input_stem, HAVE_PANDAS
from .accession import accession_base, is_local_or_pseudo_accession, looks_like_fasta_path, read_accession_list, strip_blast_suffix
from .blasttab import extract_accession_numbers_from_rows, get_topN_hits_py, header_columns, merge_blast_tsv_text, parse_blast_results_py
from .engines import (compute_parallel_params, prepare_engine_db, run_biopython_qblast_one_batch,
                      run_diamond_one_batch, run_local_blast_one_batch, run_mmseqs_one_batch,
                      run_remote_blast_one_batch)
from .fastautil import BatchInfo, detect_sequence_type, get_blast_db, get_blast_program, iter_fasta_batches
from .fields import ACCESSION_LIST_EXTENSIONS, ADD_HEADER_EXTENSIONS, BLAST_HEADERS, BLAST_REQUIRED_MODES, BLAST_STD_FIELDS, BLAST_TABLE_EXTENSIONS, NT_META_COLS, _FLOAT_FIELDS, _INT_FIELDS, build_outfmt_string, columns_from_outfmt
from .metadata import fetch_metadata_batch, load_local_meta, local_row_to_ncbi_meta, merge_rows_with_metadata
from .util import ensure_dir, is_file_complete, logger, which
from .writers import debug_local_hit_examples, rows_to_delimited, rows_to_tsv, write_results_csv

# ─────────────────────────────────────────────
# I/O discovery
# ─────────────────────────────────────────────
def find_input_files(
    input_path: str, extensions: List[str], suffix: str, mode: str
) -> List[str]:
    p = Path(input_path)
    if mode == "metadata-only":
        if p.is_file():
            return [str(p)]
        if p.is_dir():
            # Codex 修改: metadata-only 目录输入支持常见 accession-list 扩展名。
            # 理由: 原代码只查 *.txt，旧版帮助中允许 accession list，实际 .tsv/.csv 常见。
            # 原代码:
            #   return [str(x) for x in sorted(p.glob("*.txt"))]
            files: List[Path] = []
            for ext in ACCESSION_LIST_EXTENSIONS:
                files.extend(p.glob(f"*{ext}"))
            return [str(x) for x in sorted(set(files))]
        return []

    if mode == "add-header":
        if p.is_file():
            return [str(p)]
        if p.is_dir():
            found: List[Path] = []
            for ext in ADD_HEADER_EXTENSIONS:
                found.extend(x for x in p.glob(f"*{ext}") if not x.stem.endswith("_with_header"))
            return [str(x) for x in sorted(set(found))]
        return []

    if mode == "topN-only":
        # Codex 修改: topN-only 允许直接输入已有 BLAST 表，或目录中的 BLAST 表。
        # 理由: 原逻辑按 FASTA 扩展名找输入，导致 topN-only 无法从既有 BLAST 结果生成 TopN。
        if p.is_file():
            return [str(p)]
        if p.is_dir():
            files: List[Path] = []
            for ext in list(extensions) + BLAST_TABLE_EXTENSIONS:
                files.extend(p.glob(f"*{ext}"))
            return [str(x) for x in sorted(set(files))]
        return []

    # 修改 (2026-10-22): FASTA 输入接受 .gz 压缩文件。
    if p.is_file():
        ok = any(str(p).endswith(ext + suffix) or str(p).endswith(ext + suffix + ".gz") for ext in extensions)
        return [str(p)] if ok else []

    if p.is_dir():
        files: List[str] = []
        for ext in extensions:
            files.extend(glob.glob(str(p / f"*{ext}{suffix}")))
            files.extend(glob.glob(str(p / f"*{ext}{suffix}.gz")))
        return sorted(files)

    return []


def setup_output_directory(output_dir: Path) -> Dict[str, Path]:
    output_dir = ensure_dir(output_dir)
    return {
        "blast":   ensure_dir(output_dir / "blast_results"),
        "topn":    ensure_dir(output_dir / "topN_results"),
        "meta":    ensure_dir(output_dir / "metadata_results"),
        "summary": ensure_dir(output_dir / "summary_results"),
        "temp":    ensure_dir(output_dir / "temp"),
    }


def _write_completion(marker: Path, args: argparse.Namespace, seq_type: Optional[str]) -> None:
    with open(marker, "w", encoding="utf-8") as f:
        f.write(f"completed_at: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"mode: {args.mode}\n")
        if seq_type:
            f.write(f"seq_type: {seq_type}\n")
        f.write(f"blast_method: {args.blast_method}\n")
        f.write(f"threads: {args.threads}\n")
        f.write(f"parallel_jobs: {args.parallel_jobs}\n")
        f.write(f"batch_size: {args.batch_size}\n")
        f.write(f"pandas: {HAVE_PANDAS}\n")
        f.write(f"entrez_tool: {args.entrez_tool}\n")
        f.write(f"api_key_set: {bool(args.api_key)}\n")


def preflight_checks(args: argparse.Namespace, files: List[str], out_dir: Path) -> None:
    issues:   List[str] = []
    warnings: List[str] = []
    # Codex 修改: 去重 fatal issue，避免多个输入文件重复报告同一缺失程序。
    # 理由: 新增按输入类型检查 BLAST+ 可执行程序后，重复信息会降低可读性。
    issue_seen = set()

    def add_issue(msg: str) -> None:
        if msg not in issue_seen:
            issue_seen.add(msg)
            issues.append(msg)

    if not files:
        add_issue("No input FASTA / input files found.")

    for fp in files[:20]:
        p = Path(fp)
        if not p.exists():
            add_issue(f"Missing input file: {p}")
        elif p.stat().st_size == 0:
            add_issue(f"Empty input file: {p}")

    if args.mode in BLAST_REQUIRED_MODES:
        # Codex 修改: BLAST 预检按实际序列类型和 online_engine 检查 blastn/blastp。
        # 理由: 原代码只在 local 模式固定检查 nucleotide/blastn；默认 blast_remote 缺 blastn
        #       会到运行阶段才失败，并且旧 main 还返回 0。
        # 原代码:
        #   if args.blast_method == "local":
        #       program = args.blast_exe if args.blast_exe != "auto" else get_blast_program("nucleotide")
        #       ...
        if args.engine != "blast" and args.blast_method != "local":
            add_issue(f"--engine {args.engine} only works with -bm local (DIAMOND/MMseqs2 are local tools)")
        if args.blast_method == "local" and args.engine == "diamond" and not which(args.diamond_exe):
            add_issue(f"DIAMOND executable not found: {args.diamond_exe} (--diamond_exe)")
        if args.blast_method == "local" and args.engine == "mmseqs" and not which(args.mmseqs_exe):
            add_issue(f"MMseqs2 executable not found: {args.mmseqs_exe} (--mmseqs_exe)")
        for fp in files[:20]:
            seq_type = detect_sequence_type(fp)
            program = args.blast_exe if args.blast_exe != "auto" else get_blast_program(seq_type)
            needs_blast_plus = (args.blast_method == "local" and args.engine == "blast") or (
                args.blast_method == "online" and args.online_engine == "blast_remote"
            )
            if needs_blast_plus and not which(program):
                add_issue(f"BLAST+ executable not found in PATH for {seq_type} input: {program}")
            if args.blast_method == "local" and args.engine == "diamond" and seq_type == "nucleotide":
                warnings.append(f"{Path(fp).name}: nucleotide query with --engine diamond -> translated "
                                "search (blastx) against a protein database")

            if args.blast_method == "local" and args.engine == "blast":
                db = get_blast_db(seq_type, args)
                db_p = Path(str(db))
                suffixes = (
                    [".pin", ".phr", ".psq", ".00.pin", ".00.phr", ".00.psq"]
                    if seq_type == "protein"
                    else [".nin", ".nhr", ".nsq", ".00.nin", ".00.nhr", ".00.nsq"]
                )
                if db not in {"nt", "nr"} and not db_p.exists() and not any(
                    Path(str(db) + s).exists() for s in suffixes
                ):
                    warnings.append(f"Local BLAST DB path could not be verified directly: {db}")

    if args.local_meta_tsv:
        lp = Path(args.local_meta_tsv)
        if not lp.exists():
            add_issue(f"--local_meta_tsv not found: {lp}")
        elif lp.stat().st_size == 0:
            add_issue(f"--local_meta_tsv is empty: {lp}")

    if args.mode in {"blast2meta", "blast2topN2metadata", "metadata-only"}:
        # Codex 修改: 对 placeholder email 给出预警。
        # 理由: Entrez 请求需要真实 email，默认占位符容易造成合规和限流问题。
        if not args.no_ncbi_metadata and args.email in {"", "you@example.com"}:
            warnings.append(
                "NCBI metadata fetch is enabled but --email is empty or the placeholder you@example.com."
            )

    if args.use_snakemake and not which("snakemake"):
        warnings.append("snakemake not found; workflow will fall back to plain Python mode.")

    out_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "fatal_issues":  issues,
        "warnings":      warnings,
        "n_input_files": len(files),
        "input_preview": files[:20],
        "output_dir":    str(out_dir),
        "timestamp":     time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(out_dir / "preflight_report.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)

    for w in warnings:
        logger.warning(f"[PREFLIGHT] {w}")
    for e in issues:
        logger.error(f"[PREFLIGHT] {e}")

    if issues:
        raise SystemExit("Preflight checks failed. See preflight_report.json")


# ─────────────────────────────────────────────
# main per-FASTA (Fix-H: csv_output always written)
# ─────────────────────────────────────────────
def process_one_fasta(
    fasta_file: str,
    args: argparse.Namespace,
    dirs: Dict[str, Path],
) -> bool:
    fasta_path = Path(fasta_file)
    # Codex 修改: process_one_fasta 返回 bool，并识别 topN-only 的 BLAST 表输入。
    # 理由: 原代码返回 None，失败路径只 log 后 return，主程序无法产生非零退出码；
    #       topN-only 也无法直接处理已有 BLAST TSV。
    # 原签名:
    #   def process_one_fasta(... ) -> None:
    fasta_exts = [x.strip() for x in args.extension_fasta.split(",") if x.strip()]
    input_is_blast_table = args.mode == "topN-only" and not looks_like_fasta_path(fasta_path, fasta_exts)
    file_stem = (
        args.out_stem
        if args.out_stem
        else strip_blast_suffix(fasta_path.stem) if input_is_blast_table else input_stem(fasta_path)
    )

    blast_output_merged = dirs["blast"]   / f"{file_stem}_blast.tsv"
    topn_output         = dirs["topn"]    / f"{file_stem}_topN.tsv"
    metadata_output     = dirs["meta"]    / f"{file_stem}_metadata.tsv"
    summary_output      = dirs["summary"] / f"{file_stem}_summary.tsv"
    csv_output          = Path(args.output) / f"{file_stem}_results.csv"
    completion_marker   = Path(args.output) / f"{file_stem}_completed.txt"

    if args.resume and is_file_complete(csv_output) and completion_marker.exists():
        logger.info(f"[SKIP] {fasta_file} (resume: outputs exist)")
        return True

    blast_needed = args.mode in ["blast-only", "blast2topN", "blast2meta", "blast2topN2metadata"]
    topn_needed  = args.mode in ["topN-only", "blast2topN", "blast2topN2metadata"]
    meta_needed  = args.mode in ["blast2meta", "blast2topN2metadata", "metadata-only"]

    seq_type = None
    # Codex 修改: 只在真正要跑 BLAST 时检测序列类型。
    # 理由: 原代码在 topN-only 也尝试把 BLAST TSV 当 FASTA 检测，容易误导日志。
    # 原代码:
    #   if args.mode != "metadata-only":
    if blast_needed:
        seq_type = detect_sequence_type(fasta_path)
        logger.info(f"[SEQTYPE] {fasta_path.name}: {seq_type}")

    local_meta: Dict[str, Dict[str, str]] = {}
    if args.local_meta_tsv:
        local_meta, local_meta_key_used, fieldnames = load_local_meta(
            args.local_meta_tsv, args.local_meta_key
        )
        logger.info(
            f"[LOCAL_META] loaded rows={len(local_meta)} "
            f"key_requested={args.local_meta_key} key_used={local_meta_key_used}"
        )
        if fieldnames:
            logger.info(f"[LOCAL_META] columns={','.join(fieldnames)}")
        if len(local_meta) == 0:
            logger.warning(
                "[LOCAL_META] Loaded 0 rows. Check --local_meta_tsv, and whether the key "
                "matches the file. For local_nt_metadata.tsv the best key is usually "
                "accession or auto."
            )

    pseudo_prefixes = [
        x.strip() for x in (args.pseudo_acc_prefixes or "").split(",") if x.strip()
    ]

    # Fix-2/3: derive col_names once so parser and writers stay consistent
    _, col_names = build_outfmt_string(args.outfmt_blast)
    # Codex 修改: qblast 后端禁止自定义列顺序。
    # 理由: Biopython qblast 解析 XML 后固定写 14 列，不能遵守 BLAST+ outfmt 的自定义列顺序。
    fixed_columns_engine = (args.blast_method == "online" and args.online_engine == "biopython") or (
        args.blast_method == "local" and args.engine != "blast"
    )
    if blast_needed and fixed_columns_engine and col_names != BLAST_HEADERS:
        logger.error(
            "[ENGINE] biopython qblast, DIAMOND and MMseqs2 write the built-in 14 BLAST columns; "
            "custom --outfmt_blast column order is only supported with local/remote BLAST+."
        )
        return False
    if blast_needed and args.blast_method == "local" and args.engine != "blast":
        prepare_engine_db(args, seq_type, Path(args.output))

    # ── 1) BLAST ──────────────────────────────
    if blast_needed:
        if args.resume and is_file_complete(blast_output_merged):
            logger.info(f"[BLAST] merged exists, skip: {blast_output_merged.name}")
        else:
            batch_size = int(args.batch_size or 50)
            batches    = iter_fasta_batches(fasta_path, dirs["temp"] / file_stem, batch_size)
            if not batches:
                logger.warning(f"[EMPTY_FASTA] {fasta_file}")
                return False

            blast_files: List[Path] = []
            if args.blast_method == "local":
                parallel_jobs, threads_per_job = compute_parallel_params(args, mode="local", n_batches=len(batches))
                logger.info(
                    f"[LOCAL_SCHED] total_threads={args.threads} "
                    f"parallel_jobs={parallel_jobs} threads_per_job={threads_per_job}"
                )
                local_runner = {
                    "blast": run_local_blast_one_batch,
                    "diamond": run_diamond_one_batch,
                    "mmseqs": run_mmseqs_one_batch,
                }[args.engine]

                def runner(b: BatchInfo) -> Tuple[str, bool, Path, List[str]]:
                    return local_runner(
                        b, seq_type, args, dirs["blast"] / file_stem, threads_per_job, col_names
                    )
            else:
                parallel_jobs, _ = compute_parallel_params(args, mode="online")
                logger.info(
                    f"[ONLINE_SCHED] parallel_jobs={parallel_jobs} "
                    f"engine={args.online_engine} batch_size={args.batch_size} "
                    f"delay={args.online_delay_sec}s"
                )
                if args.online_engine == "blast_remote":
                    # 修改 (2026-10-10): 代理环境提示。实测 (blastn 2.17, HTTP 代理): -remote 直接 rc=255
                    # "Connection stream is in bad state"；设 CONN_HTTP_PROXY_HOST/PORT 后 20 分钟无结果；
                    # 同 query 用 biopython 引擎 <1 分钟完成。
                    if any(os.environ.get(k) for k in ("https_proxy", "HTTPS_PROXY", "http_proxy", "HTTP_PROXY")):
                        logger.warning(
                            "[ONLINE] HTTP proxy detected: blastn -remote often fails or hangs behind a proxy; "
                            "if it does, use --online_engine biopython"
                        )
                    def runner(b: BatchInfo) -> Tuple[str, bool, Path, List[str]]:
                        return run_remote_blast_one_batch(
                            b, seq_type, args, dirs["blast"] / file_stem, col_names
                        )
                else:
                    def runner(b: BatchInfo) -> Tuple[str, bool, Path, List[str]]:
                        return run_biopython_qblast_one_batch(
                            b, seq_type, args, dirs["blast"] / file_stem, col_names
                        )

            # 修改 (2026-10-22): Ctrl-C 时取消尚未开始的批次，原代码 with 块退出会等所有批次跑完。
            with cf.ThreadPoolExecutor(
                max_workers=min(parallel_jobs, len(batches))
            ) as ex:
                futs = [ex.submit(runner, b) for b in batches]
                try:
                    for fut in cf.as_completed(futs):
                        bid, ok, outp, _ = fut.result()
                        if ok:
                            blast_files.append(outp)
                        else:
                            logger.error(f"[BLAST_FAIL] batch={bid}")
                except KeyboardInterrupt:
                    for f in futs:
                        f.cancel()
                    raise

            # 修改 (2026-10-10): 任一批次失败即判本文件失败，不合并、不写完成标记。
            # 理由: 原代码只 log 后合并成功批次并返回成功；缺失 query 无提示，且合并文件存在后
            #       --resume 不会重跑失败批次。成功批次结果保留，--resume 只重跑失败批次。
            # 原代码:
            #   if not blast_files:
            #       logger.error(...); return False
            n_failed = len(batches) - len(blast_files)
            if n_failed:
                logger.error(
                    f"[BLAST] {n_failed}/{len(batches)} batches failed for {fasta_file}; "
                    "rerun with --resume to retry only the failed batches"
                )
                return False

            merge_blast_tsv_text(
                blast_files, blast_output_merged, None if args.no_blast_header else col_names
            )
            logger.info(f"[BLAST] merged -> {blast_output_merged}")

        # Fix-H.1: blast-only mode must write csv_output
        if args.mode == "blast-only":
            merged_rows = parse_blast_results_py(blast_output_merged, col_names)
            write_results_csv(csv_output, merged_rows, args, default_headers=col_names)
            logger.info(f"[CSV] wrote -> {csv_output}")
            _write_completion(completion_marker, args, seq_type)
            return True

    merged_rows: List[Dict[str, object]] = []
    if blast_needed and blast_output_merged.exists():
        merged_rows = parse_blast_results_py(blast_output_merged, col_names)
    elif args.mode == "topN-only":
        # Codex 修改: topN-only 从已有 BLAST 表加载 merged_rows。
        # 理由: 原代码中 blast_needed=False，merged_rows 永远为空，导致静默输出空 TopN。
        # 原代码:
        #   merged_rows = []
        #   if blast_needed and blast_output_merged.exists(): ...
        blast_source = fasta_path if input_is_blast_table else blast_output_merged
        if not is_file_complete(blast_source):
            logger.error(
                f"[TOPN] BLAST results not found for topN-only mode: {blast_source}. "
                "Use a BLAST TSV as -i or keep the prior output under blast_results/{stem}_blast.tsv."
            )
            return False
        merged_rows = parse_blast_results_py(blast_source, col_names)
        logger.info(f"[TOPN] loaded BLAST table -> {blast_source}")

    if local_meta and merged_rows:
        debug_local_hit_examples(merged_rows, local_meta, pseudo_prefixes, n=8)

    # ── 2) TopN ──────────────────────────────
    top_rows: Optional[List[Dict[str, object]]] = None
    if topn_needed:
        if args.resume and is_file_complete(topn_output):
            logger.info(f"[TOPN] exists, skip: {topn_output.name}")
            top_rows = parse_blast_results_py(topn_output, col_names)
        else:
            if not merged_rows:
                logger.warning(f"[TOPN] no hits for {fasta_file}")
                top_rows = []
            else:
                top_rows = get_topN_hits_py(
                    merged_rows, int(args.topN), args.sort_by,
                    min_pident=float(args.min_pident), min_qcovs=float(args.min_qcovs),
                )
            rows_to_tsv(topn_output, top_rows or [], default_headers=col_names)
            logger.info(f"[TOPN] wrote -> {topn_output}")

        # Fix-H.2: topN-only mode must write csv_output directly
        if args.mode == "topN-only":
            write_results_csv(csv_output, top_rows or [], args, default_headers=col_names)
            logger.info(f"[CSV] wrote -> {csv_output}")
            _write_completion(completion_marker, args, seq_type)
            return True

        if args.mode == "blast2topN":
            # Codex 修改: blast2topN 也写最终 *_results.csv。
            # 理由: 原代码只写 topN_results/*.tsv，Snakemake rule all 和 --resume 期待 CSV 哨兵。
            write_results_csv(csv_output, top_rows or [], args, default_headers=col_names)
            logger.info(f"[CSV] wrote -> {csv_output}")
            _write_completion(completion_marker, args, seq_type)
            return True

    # ── 3) Metadata fetch + merge ─────────────
    if meta_needed:
        metadata_dict: Dict[str, dict] = {}
        # Codex 修改: metadata-only 读取 accession list；其他模式仍从 BLAST/TopN 行提取。
        # 理由: 原代码在 metadata-only 下直接 accessions=[]，该模式实际空跑。
        # 原代码:
        #   if args.mode == "metadata-only":
        #       accessions: List[str] = []
        if args.mode == "metadata-only":
            accessions = read_accession_list(fasta_path)
            logger.info(f"[META_ONLY] loaded accessions={len(accessions)} from {fasta_path}")
        elif args.mode == "blast2meta":
            accessions = extract_accession_numbers_from_rows(merged_rows)
        else:
            accessions = extract_accession_numbers_from_rows(top_rows or [])

        if args.mode == "metadata-only" and local_meta:
            # Codex 修改: metadata-only + local_meta_tsv 时可直接输出本地 metadata。
            # 理由: 旧逻辑只有 BLAST summary merge 才能使用 local_meta，metadata-only 无法替代本地库用法。
            for acc in accessions:
                acc_full = (acc or "").strip()
                acc_base = accession_base(acc_full)
                lrow = local_meta.get(acc_full) or local_meta.get(acc_base)
                if not lrow:
                    continue
                if args.local_meta_mode == "map_to_ncbi":
                    lm = local_row_to_ncbi_meta(lrow)
                    key = lm.get("accession") or acc_base
                    metadata_dict[key] = lm
                else:
                    metadata_dict[acc_base] = {
                        "accession": acc_base,
                        **{args.local_meta_prefix + k: v for k, v in lrow.items()},
                    }

        if not args.no_ncbi_metadata:
            accessions_real = [
                a for a in accessions
                if not is_local_or_pseudo_accession(a, local_meta, pseudo_prefixes)
            ]
            logger.info(
                f"[META] accessions(total)={len(accessions)} "
                f"real_for_entrez={len(accessions_real)} "
                f"local_or_pseudo_skipped={len(accessions) - len(accessions_real)}"
            )

            if args.resume and is_file_complete(metadata_output):
                logger.info(f"[META] exists, load: {metadata_output.name}")
                with open(metadata_output, "r", encoding="utf-8", errors="ignore") as f:
                    reader = csv.DictReader(f, delimiter="\t")
                    for r in reader:
                        a = str(r.get("accession", "")).split(".")[0]
                        if a:
                            metadata_dict[a] = dict(r)

            missing = [a for a in accessions_real if a not in metadata_dict]
            if missing:
                preferred_db = "protein" if seq_type == "protein" else "nucleotide"
                metadata_new = fetch_metadata_batch(
                    missing, args, batch_size=20, preferred_db=preferred_db
                )
                metadata_dict.update(metadata_new)

            if metadata_dict:
                rows_to_tsv(metadata_output, list(metadata_dict.values()), default_headers=NT_META_COLS)
                logger.info(f"[META] wrote -> {metadata_output}")
        else:
            logger.info("[META] --no_ncbi_metadata ON, skip Entrez efetch.")
            if metadata_dict:
                rows_to_tsv(metadata_output, list(metadata_dict.values()), default_headers=NT_META_COLS)
                logger.info(f"[META] wrote local metadata -> {metadata_output}")

        # Fix-H.3: metadata-only mode must write csv_output directly
        if args.mode == "metadata-only":
            write_results_csv(csv_output, list(metadata_dict.values()), args, default_headers=NT_META_COLS)
            logger.info(f"[CSV] wrote -> {csv_output}")
            _write_completion(completion_marker, args, seq_type)
            return True

        base_rows = merged_rows if args.mode == "blast2meta" else (top_rows or [])
        merged    = merge_rows_with_metadata(base_rows, metadata_dict, local_meta, args)

        rows_to_tsv(summary_output, merged)
        # Fix-H.4: blast2*/no-pandas fallback must write csv_output directly
        write_results_csv(csv_output, merged, args)

        logger.info(f"[SUMMARY] wrote -> {summary_output}")
        logger.info(f"[CSV] wrote -> {csv_output}")

    _write_completion(completion_marker, args, seq_type)
    return True



# ─────────────────────────────────────────────
# add-header mode (2026-10-11)
# ─────────────────────────────────────────────

def _check_column_types(rows: List[List[str]], cols: List[str]) -> Optional[str]:
    """Return a message if values do not fit the column types, else None."""
    for row in rows[:200]:
        for code, v in zip(cols, row):
            v = v.strip()
            if code in _INT_FIELDS or code in _FLOAT_FIELDS:
                try:
                    x = float(v)
                except ValueError:
                    return f"column '{code}' has non-numeric value {v!r}"
                if code in _INT_FIELDS and not x.is_integer():
                    return f"column '{code}' should be an integer, got {v!r}"
                if code in {"pident", "ppos", "qcovs", "qcovhsp", "qcovus"} and not 0 <= x <= 100:
                    return f"column '{code}' should be 0-100, got {v!r}"
    return None


def add_header_to_table(table: Path, args: argparse.Namespace, out_dir: Path) -> bool:
    """Write <stem>_with_header.csv with a header row for an existing BLAST tabular file.

    Column source, in order: outfmt 7 '# Fields:' line or an existing header row in the file;
    --outfmt_blast with an explicit column list (e.g. "6 std staxids"); otherwise guessed from
    the column count (12 = BLAST default 'std', 14 = this script's default) and type-checked.
    """
    cols: Optional[List[str]] = None
    source = ""
    data_lines: List[str] = []
    with open(table, "r", encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            fields = header_columns(line)
            if fields:
                if cols is None:
                    cols, source = fields, "header line in file"
                continue
            if line.strip() and not line.startswith("#"):
                data_lines.append(line.rstrip("\r\n"))
    if not data_lines:
        logger.warning(f"[ADD_HEADER] no data rows: {table}")
    delim = "\t" if not data_lines or "\t" in data_lines[0] else ","
    rows = [r for r in csv.reader(data_lines, delimiter=delim)] if delim == "," else [x.split("\t") for x in data_lines]

    explicit = (args.outfmt_blast or "").strip() not in {"", "6", "7"}
    if cols is None and explicit:
        cols, source = columns_from_outfmt(args.outfmt_blast), f"--outfmt_blast '{args.outfmt_blast}'"
        if cols is None:
            logger.error(f"[ADD_HEADER] --outfmt_blast needs a column list, e.g. '6 std staxids': {args.outfmt_blast}")
            return False
    guessed = False
    if cols is None:
        n = len(rows[0]) if rows else 0
        guess = {12: BLAST_STD_FIELDS, 14: BLAST_HEADERS}.get(n)
        if guess is None:
            logger.error(
                f"[ADD_HEADER] {table.name}: {n} columns and no header/'# Fields:' line; cannot know the "
                "columns. Pass the -outfmt you used, e.g. --outfmt_blast \"6 qseqid sseqid pident evalue\""
            )
            return False
        cols, guessed = list(guess), True
        source = "BLAST default 'std' (12 columns)" if n == 12 else "blast2_metadata default (14 columns)"

    width = len(cols)
    fixed: List[List[str]] = []
    for r in rows:
        if len(r) == width:
            fixed.append(r)
        elif len(r) > width and "stitle" in cols:
            i = cols.index("stitle")
            tail = width - i - 1
            fixed.append(r[:i] + [delim.join(r[i:len(r) - tail])] + (r[len(r) - tail:] if tail else []))
        else:
            logger.error(f"[ADD_HEADER] {table.name}: row has {len(r)} columns, expected {width} ({source}): {r[:4]}")
            return False
    problem = _check_column_types(fixed, cols)
    if problem:
        hint = " The guess is wrong; pass the -outfmt you used via --outfmt_blast." if guessed else ""
        logger.error(f"[ADD_HEADER] {table.name}: {problem} ({source}).{hint}")
        return False
    if guessed:
        logger.warning(f"[ADD_HEADER] {table.name}: no header in file; assumed {source}. "
                       "If you ran BLAST with a custom -outfmt, pass it via --outfmt_blast.")

    out = out_dir / f"{strip_blast_suffix(table.stem)}_with_header.csv"
    dict_rows: List[Dict[str, object]] = [dict(zip(cols, r)) for r in fixed]
    rows_to_delimited(out, dict_rows, ",", default_headers=cols,
                      header_style=args.header_style, encoding="utf-8-sig")
    logger.info(f"[ADD_HEADER] {table.name}: {len(fixed)} rows, columns from {source} -> {out}")
    return True
