#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MitoFinder 整合处理工具 v2.6
==============================
三种运行模式（subcommand 风格）：

  only_assembly      仅对 reads 进行批量组装（megahit / metaspades）
                     → 输出 contig fasta，不做线粒体注释

  assembly2summary   reads → MitoFinder（组装+注释）→ 提取汇总  全流程

  annotation2summary 已有 contig fasta → MitoFinder 注释 → 提取汇总

辅助子命令：
  snakemake          生成 Snakemake 工作流（支持 HPC 集群投递）
  analyze            智能分析输入目录结构，给出参数建议

所有路径均支持绝对路径和相对路径，入口统一 resolve 为绝对路径。
Singularity 运行时自动挂载所需目录，保证容器内路径一致。

v2.6 Bug 修复：
  BUG-K  only_assembly 断点重投脏目录必崩溃（Dirty Output Dir Crash）
         megahit / SPAdes 发现 --output/-o 目录已存在时立即退出报错：
           "MEGAHIT Error: output directory already exists"
           "Error: the output directory already exists"
         v2.5 虽删除了 ensure_dir(sample_out) 预创建，但重新运行时
         残留目录依然存在，触发同样崩溃。
         修复：组装器启动前主动 shutil.rmtree 清理已存在的 sample_out；
         并记录警告日志，让用户感知清理行为。

  BUG-L  Singularity 回退机制误伤无镜像用户（Image Not Found Crash）
         原逻辑：mitofinder_script_path 为空时无条件调用 Singularity，
         但若用户机器上 Singularity 已安装而默认 SIF
         （/home/deyuan/Softs/mitofinder_v1.4.1.sif）不存在，
         则命令以 "image not found" 崩溃，且这个错误极难定位。
         影响范围：_run_assembly_only_task（only_assembly 模式）
                   _run_single_task（assembly2summary / annotation2summary）
         修复：新增三级判定逻辑
           1. 用户显式指定 --singularity_sif → 必定使用容器
           2. singularity 命令存在 且 默认 SIF 文件实际存在 → 顺手用容器
           3. 否则 → 直接调用宿主机命令（inner_cmd），并输出 INFO 提示

v2.5 Bug 修复（继承）：
  BUG-F  only_assembly：删除 ensure_dir(sample_out) 预创建
  BUG-G  _build_mf_args_assembly：删除无效的 -l 参数
  BUG-H  collect_assembly_tasks TSV 相对路径解析
  BUG-I  record_script 使用 shlex.join 保留参数引号
  BUG-J  gather_genes 重名序列 ID 去重

继承 v2.2–v2.4 全部修复（FIX-1~6, BUG-A~E）。
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import textwrap
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Generator, List, Optional, Set, Tuple

__version__ = "2.6.0"

# ── 可选依赖 ──────────────────────────────────────────────────────────────────
try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("[ERROR] Biopython 未安装，请执行: pip install biopython", file=sys.stderr)
    sys.exit(1)

# ══════════════════════════════════════════════════════════════════════════════
# 常量
# ══════════════════════════════════════════════════════════════════════════════

_DEFAULT_SIF = "/home/deyuan/Softs/mitofinder_v1.4.1.sif"
_OFFICIAL_MITOFINDER_URL = "https://github.com/RemiAllio/MitoFinder"

# ══════════════════════════════════════════════════════════════════════════════
# 日志  [FIX-1] 支持动态追加 FileHandler
# ══════════════════════════════════════════════════════════════════════════════

def setup_logging(log_file: Optional[Path] = None, verbose: bool = False) -> logging.Logger:
    logger = logging.getLogger("mitofinder")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    fmt = logging.Formatter("[%(asctime)s %(levelname)s] %(message)s", datefmt="%H:%M:%S")

    has_stream = any(
        isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
        for h in logger.handlers
    )
    if not has_stream:
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(fmt)
        logger.addHandler(ch)

    if log_file:
        target = str(log_file.resolve())
        has_file = any(
            isinstance(h, logging.FileHandler) and h.baseFilename == target
            for h in logger.handlers
        )
        if not has_file:
            log_file.parent.mkdir(parents=True, exist_ok=True)
            fh = logging.FileHandler(log_file, encoding="utf-8")
            fh.setFormatter(fmt)
            logger.addHandler(fh)

    return logger


log: logging.Logger = logging.getLogger("mitofinder")

# ══════════════════════════════════════════════════════════════════════════════
# 数据类
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class AssemblyTask:
    """单个样本的组装任务，所有路径均为绝对路径。"""
    seq_id:   str
    pe1:      Optional[Path] = None   # forward reads
    pe2:      Optional[Path] = None   # reverse reads
    se:       Optional[Path] = None   # single-end reads
    assembly: Optional[Path] = None   # 已有组装 fasta（跳过独立组装步骤）


@dataclass
class AnnotateTask:
    """单个 MitoFinder 注释任务（已有 contig fasta），所有路径为绝对路径。"""
    fasta_path: Path
    seq_id:     str


@dataclass
class DirAnalysis:
    structure:       str       = "unknown"
    toplevel_fastas: List[str] = field(default_factory=list)
    subdirs:         List[str] = field(default_factory=list)
    spades_dirs:     List[str] = field(default_factory=list)
    recommended:     dict      = field(default_factory=dict)

# ══════════════════════════════════════════════════════════════════════════════
# 公共工具函数
# ══════════════════════════════════════════════════════════════════════════════

_finished_lock = threading.Lock()


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def safe_resolve(path_str: str, must_exist: bool = False) -> Path:
    """将字符串路径（绝对或相对）解析为绝对 Path。"""
    p = Path(path_str).expanduser().resolve()
    if must_exist and not p.exists():
        raise FileNotFoundError(f"路径不存在: {p}")
    return p


def load_finished(finished_file: Path) -> Set[str]:
    if not finished_file.exists():
        return set()
    return {
        line.strip()
        for line in finished_file.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }


def record_finished(finished_file: Path, seq_id: str) -> None:
    with _finished_lock:
        with finished_file.open("a", encoding="utf-8") as f:
            f.write(seq_id + "\n")


def record_script(script_file: Path, cmd_line: str) -> None:
    with _finished_lock:
        with script_file.open("a", encoding="utf-8") as f:
            f.write(cmd_line + "\n\n")


def match_pattern(name: str, pattern: Optional[str]) -> bool:
    """fnmatch 模式匹配；pattern 为 None 时返回 True（不过滤）。"""
    if not pattern:
        return True
    return fnmatch.fnmatch(name, pattern)


def _section(title: str) -> None:
    bar = "─" * 68
    log.info("\n%s\n  %s\n%s", bar, title, bar)


def iter_final_results_dirs(
    base_dir: Path,
    suffix_subdir: Optional[str] = None,
) -> Generator[Tuple[Path, Path, str], None, None]:
    """
    遍历 base_dir，找出所有含 *_Final_Results 子目录的样本目录。
    yield (sample_dir, final_results_dir, sample_name)
    """
    if not base_dir.is_dir():
        log.warning("提取输入目录不存在: %s", base_dir)
        return
    for sample_dir in sorted(base_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        sample_name = sample_dir.name
        if suffix_subdir and not match_pattern(sample_name, suffix_subdir):
            continue
        for sub in sample_dir.iterdir():
            if sub.is_dir() and sub.name.endswith("_Final_Results"):
                yield sample_dir, sub, sample_name
                break

# ══════════════════════════════════════════════════════════════════════════════
# Singularity 路径绑定
# ══════════════════════════════════════════════════════════════════════════════

def _build_singularity_binds(paths: List[Path]) -> List[str]:
    """
    生成去重后的 Singularity --bind 参数列表（宿主机 1:1 映射）。
    [FIX-6] 检测保留字符 ':' ','。
    """
    dirs: Set[Path] = set()
    for p in paths:
        if ":" in str(p) or "," in str(p):
            log.warning("路径含 Singularity 保留字符(':' 或 ',')，可能挂载失败: %s", p)
        if not p.exists():
            log.warning("绑定路径不存在，跳过: %s", p)
            continue
        dirs.add(p if p.is_dir() else p.parent)

    sorted_dirs = sorted(dirs)
    unique: List[Path] = []
    for d in sorted_dirs:
        if not any(d != u and str(d).startswith(str(u) + os.sep) for u in unique):
            unique.append(d)
    return [f"{d}:{d}" for d in unique]


def build_singularity_cmd(
    sif_path: str,
    inner_cmd: List[str],
    bind_paths: List[Path],
    workdir: Optional[Path] = None,
) -> List[str]:
    binds = _build_singularity_binds(bind_paths)
    cmd = ["singularity", "run"]
    for b in binds:
        cmd += ["--bind", b]
    if workdir:
        cmd += ["--pwd", str(workdir)]
    cmd.append(sif_path)
    cmd.extend(inner_cmd)
    return cmd

# ══════════════════════════════════════════════════════════════════════════════
# MitoFinder 参数构建
# ══════════════════════════════════════════════════════════════════════════════

def _build_mf_args_annotate(
    seq_id: str,
    assembly: Path,
    genbank_ref: Path,
    args: argparse.Namespace,
) -> List[str]:
    """
    构建 MitoFinder 注释模式参数（-a 已有组装）。

    [OFFICIAL-MITOFINDER] 官方 README 说明 -a/--assembly 接收已有
    assembly.fasta，且该 fasta 可以包含 one or several contigs。
    因此本脚本将同一样本的多条 contig 合并到一个样本 fasta 后再传给 -a。
    """
    mf = [
        "-j", seq_id,
        "-a", str(assembly),
        "-r", str(genbank_ref),
        "-t", args.tool_param,
        "-o", str(args.genetic_code),
        "-p", str(args.threads),
        "-m", str(args.memory),
    ]
    _append_common_mf_flags(mf, args)
    return mf


def _build_mf_args_assembly(
    task: AssemblyTask,
    genbank_ref: Path,
    args: argparse.Namespace,
) -> List[str]:
    """
    构建 MitoFinder 组装+注释模式参数（PE / SE reads 或已有 assembly）。

    [BUG-A] MitoFinder 参数互斥规则：
      - reads 输入（-1/-2 或 -s）：必须同时指定组装器 Flag（--megahit 等）
      - 已有 assembly 输入（-a）：绝对不能附加组装器 Flag，否则参数冲突崩溃

    [BUG-G] MitoFinder 不存在 -l 参数，传入会立即 "Unknown option: -l" 崩溃。
      contig 最小长度过滤通过 --min-contig-size 传递，已在 _append_common_mf_flags
      中正确处理，这里不需要也不能再传 -l。
    """
    mf = [
        "-j", task.seq_id,
        "-r", str(genbank_ref),
        "-t", args.tool_param,
        "-o", str(args.genetic_code),
        "-p", str(args.threads),
        "-m", str(args.memory),
    ]

    assembler = getattr(args, "assembler", "megahit")

    if task.pe1 and task.pe2:
        # PE reads → 需要组装器 Flag
        mf += ["-1", str(task.pe1), "-2", str(task.pe2), f"--{assembler}"]
    elif task.se:
        # SE reads → 需要组装器 Flag
        mf += ["-s", str(task.se), f"--{assembler}"]
    elif task.assembly:
        # 已有 assembly → 仅注释，不附加组装器 Flag
        mf += ["-a", str(task.assembly)]

    _append_common_mf_flags(mf, args)
    return mf


def _append_common_mf_flags(mf: List[str], args: argparse.Namespace) -> None:
    """附加所有模式共用的可选 MitoFinder 标志。"""
    if getattr(args, "adjust_direction", True):
        mf.append("--adjust-direction")
    if getattr(args, "override", False):
        mf.append("--override")
    if getattr(args, "ignore_non_std", False):
        mf.append("--ignore")
    if getattr(args, "new_genes", False):
        mf.append("--new-genes")
    if getattr(args, "allow_intron", False):
        mf.append("--allow-intron")
    if getattr(args, "out_gb", False):
        mf.append("--out-gb")
    if getattr(args, "cds_merge", False):
        mf.append("--cds-merge")
    for flag, opt in [
        ("blast_eval",           "-e"),
        ("nwalk",                "-n"),
        ("min_contig_size",      "--min-contig-size"),
        ("max_contig_size",      "--max-contig-size"),
        ("circular_offset",      "--circular-offset"),
        ("circular_size",        "--circular-size"),
        ("blast_identity_nucl",  "--blast-identity-nucl"),
        ("max_contig",           "--max-contig"),
        ("intron_size",          "--intron-size"),
        ("rename_contig",        "--rename-contig"),
        ("config_file",          "-c"),
    ]:
        val = getattr(args, flag, None)
        if val is not None:
            mf += [opt, str(val)]

# ══════════════════════════════════════════════════════════════════════════════
# 通用任务执行器
# ══════════════════════════════════════════════════════════════════════════════

def _run_single_task(
    seq_id: str,
    mf_args: List[str],
    bind_paths: List[Path],
    output_dir: Path,
    args: argparse.Namespace,
    script_file: Path,
    finished_file: Path,
    finished_set: Set[str],
) -> bool:
    """
    线程安全地执行单个 MitoFinder 任务。
    返回 True 表示本次成功（或已跳过）。
    """
    with _finished_lock:
        if seq_id in finished_set:
            log.info("[SKIP] %s 已完成", seq_id)
            return True

    ensure_dir(output_dir / seq_id)
    sif_path = str(safe_resolve(args.singularity_sif)) if args.singularity_sif else _DEFAULT_SIF

    if args.mitofinder_script_path:
        cmd = [args.mitofinder_script_path] + mf_args
    else:
        # [BUG-L] 三级判定：避免默认 SIF 不存在时的 "image not found" 崩溃
        use_container = False
        if args.singularity_sif:
            use_container = True   # 用户显式指定，必须使用
        elif shutil.which("singularity") and os.path.exists(sif_path):
            use_container = True   # 命令存在且镜像文件实际存在
        else:
            log.info("未检测到可用的 Singularity 环境，直接调用宿主机 MitoFinder: mitofinder")

        if use_container:
            cmd = build_singularity_cmd(
                sif_path=sif_path,
                inner_cmd=mf_args,
                bind_paths=bind_paths + [output_dir],
                workdir=output_dir,
            )
        else:
            cmd = ["mitofinder"] + mf_args

    record_script(script_file, shlex.join(cmd))

    if args.dry_run:
        log.info("[DRY-RUN] %s\n  %s", seq_id, shlex.join(cmd))
        return True

    log.info("[START] %s", seq_id)
    try:
        result = subprocess.run(
            cmd,
            cwd=str(output_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=getattr(args, "task_timeout", None),
        )
        if result.returncode == 0:
            log.info("[SUCCESS] %s", seq_id)
            record_finished(finished_file, seq_id)
            with _finished_lock:
                finished_set.add(seq_id)
            return True
        else:
            log.error(
                "[FAIL] %s (code=%d)\nSTDOUT:\n%s\nSTDERR:\n%s",
                seq_id, result.returncode,
                result.stdout.strip()[-1500:],
                result.stderr.strip()[-1500:],
            )
    except subprocess.TimeoutExpired:
        log.error("[TIMEOUT] %s 超过 %ss 已中止", seq_id, args.task_timeout)
    except FileNotFoundError as exc:
        log.error("[NOT-FOUND] %s: %s", seq_id, exc)
        _log_mitofinder_runtime_hint()
    except Exception as exc:
        log.exception("[EXCEPTION] %s: %s", seq_id, exc)
    return False


def _mitofinder_runtime_ready(args: argparse.Namespace) -> bool:
    """
    [ADD-RUNTIME-HINT] 运行 MitoFinder 前做一次轻量预检查。
    目的：避免未安装 mitofinder 时并发任务全部启动后重复报 [NOT-FOUND]。
    """
    if args.mitofinder_script_path:
        p = Path(args.mitofinder_script_path).expanduser()
        return p.exists() or shutil.which(args.mitofinder_script_path) is not None
    if args.singularity_sif:
        return Path(args.singularity_sif).expanduser().exists()
    if shutil.which("mitofinder"):
        return True
    return shutil.which("singularity") is not None and os.path.exists(_DEFAULT_SIF)


def _log_mitofinder_runtime_hint() -> None:
    """
    [ADD-RUNTIME-HINT] 给出找不到 MitoFinder 时的明确参数用法。
    """
    log.error(
        "未找到 MitoFinder 可执行环境，任务池未启动，避免批量样本逐个失败。\n"
        "  方式1: 安装 MitoFinder 并确保命令在 PATH 中: which mitofinder\n"
        "  方式2: 指定本地可执行文件: --mitofinder_script_path /path/to/mitofinder\n"
        "  方式3: 指定 Singularity 镜像: --singularity_sif /path/to/mitofinder_v1.4.2.sif\n"
        "        官方 Singularity 镜像可用以下命令获取:\n"
        "        singularity pull --arch amd64 library://remiallio/default/mitofinder:v1.4.2\n"
        "  示例: python3 mitofinder_pipeline_v2_7.py annotation2summary "
        "-i rename_summary_all.fasta -o ann -r ref.gb "
        "--mitofinder_script_path /home/user/MitoFinder/mitofinder\n"
        f"  官方说明: {_OFFICIAL_MITOFINDER_URL}"
    )


def _run_pool(
    tasks_iter,           # iterable of (seq_id, mf_args, bind_paths)
    output_dir: Path,
    args: argparse.Namespace,
    script_file: Path,
    finished_file: Path,
    finished_set: Set[str],
) -> Tuple[int, int]:
    """并发执行一批任务，返回 (成功数, 失败数)。"""
    total_threads = args.max_workers * args.threads
    log.info("并发配置: max_workers=%d × threads/task=%d → 总线程 %d",
             args.max_workers, args.threads, total_threads)

    ok = fail = 0
    with ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        futures = {
            pool.submit(
                _run_single_task,
                seq_id, mf_args, bind_paths,
                output_dir, args, script_file, finished_file, finished_set,
            ): seq_id
            for seq_id, mf_args, bind_paths in tasks_iter
        }
        for fut in as_completed(futures):
            exc = fut.exception()
            if exc:
                log.error("任务 %s 未捕获异常: %s", futures[fut], exc)
                fail += 1
            elif fut.result():
                ok += 1
            else:
                fail += 1

    log.info("执行完毕: 成功=%d  失败=%d", ok, fail)
    return ok, fail

# ══════════════════════════════════════════════════════════════════════════════
# ── 独立组装（only_assembly）──  使用外部组装器，不调用 MitoFinder
# ══════════════════════════════════════════════════════════════════════════════

_ASSEMBLER_CMDS = {
    "megahit":    "megahit",
    "metaspades": "spades.py",
    "idba":       "idba_ud",
}


def _build_assembler_cmd(
    task: AssemblyTask,
    out_dir: Path,
    args: argparse.Namespace,
) -> List[str]:
    """
    构建独立组装命令（megahit / metaspades）。
    out_dir 已是样本专属输出目录。

    [BUG-C] idba_ud 需要预先将 FASTQ 转换为交织 FASTA，
    当前版本未实现此转换，直接抛出 NotImplementedError。
    请改用 megahit 或 metaspades。
    """
    assembler = args.assembler
    threads   = str(args.threads)
    memory_gb = str(args.memory)

    if assembler == "megahit":
        cmd = [
            "megahit",
            "-t", threads,
            "--memory", str(int(memory_gb) * 1_000_000_000),
            "-o", str(out_dir),
        ]
        if getattr(args, "min_contig_size", None):
            cmd += ["--min-contig-len", str(args.min_contig_size)]
        if task.pe1 and task.pe2:
            cmd += ["-1", str(task.pe1), "-2", str(task.pe2)]
        elif task.se:
            cmd += ["-r", str(task.se)]

    elif assembler == "metaspades":
        cmd = [
            "spades.py", "--meta",
            "-t", threads,
            "-m", memory_gb,
            "-o", str(out_dir),
        ]
        if task.pe1 and task.pe2:
            cmd += ["-1", str(task.pe1), "-2", str(task.pe2)]
        elif task.se:
            cmd += ["-s", str(task.se)]

    elif assembler == "idba":
        # [BUG-C] idba_ud 只接受单个交织 FASTA，需要预先用 fq2fa 转换，
        # 当前版本未实现此步骤，明确报错提示用户。
        raise NotImplementedError(
            "only_assembly --assembler idba：\n"
            "  idba_ud 需要预先将 FASTQ 合并转换为交织 FASTA（fq2fa --merge），\n"
            "  当前版本暂未实现此自动转换。\n"
            "  请改用 --assembler megahit 或 --assembler metaspades，\n"
            "  或手动转换后使用 --reads_tsv 的 assembly 列输入已有 FASTA。"
        )
    else:
        raise ValueError(f"未知组装器: {assembler}")

    return cmd


def _run_assembly_only_task(
    task: AssemblyTask,
    output_dir: Path,
    args: argparse.Namespace,
    script_file: Path,
    finished_file: Path,
    finished_set: Set[str],
) -> bool:
    """
    执行单样本独立组装（不调用 MitoFinder 注释）。

    [BUG-F] 绝对不能预先 ensure_dir(sample_out)：
      megahit 和 SPAdes 都强制要求 -o 指定的目录在运行时不存在，
      若目录已存在则直接报错退出（megahit: "output directory already exists"）。
      修复：只确保父目录 output_dir 存在，由组装器自己创建 sample_out。
      Singularity 挂载父目录 output_dir，不挂载尚未创建的 sample_out。

    [BUG-C] 组装命令套入 Singularity 容器执行，保证环境一致。
    """
    with _finished_lock:
        if task.seq_id in finished_set:
            log.info("[SKIP] %s 已完成", task.seq_id)
            return True

    # [BUG-F] 只创建父目录，sample_out 由组装器自行创建
    ensure_dir(output_dir)
    sample_out = output_dir / task.seq_id  # 此时不创建

    # [BUG-K] 断点重投保护：如果上次失败的残留目录仍然存在，
    # megahit / SPAdes 都会以 "output directory already exists" 立即崩溃。
    # 在组装器启动前主动清理，让重投可以正常进行。
    if sample_out.exists():
        log.warning("[CLEANUP] 检测到已存在的组装目录，清理以允许重新运行: %s", sample_out)
        shutil.rmtree(sample_out, ignore_errors=True)

    # [BUG-C] idba 在 _build_assembler_cmd 中会抛出 NotImplementedError
    try:
        inner_cmd = _build_assembler_cmd(task, sample_out, args)
    except NotImplementedError as exc:
        log.error("[SKIP] %s: %s", task.seq_id, exc)
        return False

    sif_path = str(safe_resolve(args.singularity_sif)) if args.singularity_sif else _DEFAULT_SIF

    if args.mitofinder_script_path:
        # 用户显式指定了本地可执行路径，直接调用宿主机组装器
        cmd = inner_cmd
        log.debug("使用本地组装器: %s", inner_cmd[0])
    else:
        # [BUG-L] 三级判定：
        #   1. 用户显式指定 --singularity_sif → 必定使用容器
        #   2. singularity 命令存在 且 默认 SIF 文件实际存在 → 使用容器
        #   3. 否则 → 直接调用宿主机命令，避免 "image not found" 崩溃
        use_container = False
        if args.singularity_sif:
            use_container = True   # 用户显式指定，必须使用
        elif shutil.which("singularity") and os.path.exists(sif_path):
            use_container = True   # 命令存在且镜像文件实际存在
        else:
            log.info("未检测到可用的 Singularity 环境，直接调用宿主机组装器: %s", inner_cmd[0])

        if use_container:
            # [BUG-F+L] 挂载父目录 output_dir（sample_out 尚未存在）
            bind_paths = [p for p in [task.pe1, task.pe2, task.se, output_dir] if p]
            cmd = build_singularity_cmd(
                sif_path=sif_path,
                inner_cmd=inner_cmd,
                bind_paths=bind_paths,
                workdir=output_dir,
            )
        else:
            cmd = inner_cmd

    record_script(script_file, shlex.join(cmd))

    if args.dry_run:
        log.info("[DRY-RUN] %s\n  %s", task.seq_id, shlex.join(cmd))
        return True

    log.info("[START-ASSEMBLY] %s  assembler=%s", task.seq_id, args.assembler)
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=getattr(args, "task_timeout", None),
        )
        if result.returncode == 0:
            log.info("[SUCCESS-ASSEMBLY] %s", task.seq_id)
            record_finished(finished_file, task.seq_id)
            with _finished_lock:
                finished_set.add(task.seq_id)
            return True
        else:
            log.error("[FAIL-ASSEMBLY] %s (code=%d)\nSTDOUT:\n%s\nSTDERR:\n%s",
                      task.seq_id, result.returncode,
                      result.stdout.strip()[-800:],
                      result.stderr.strip()[-1200:])
    except subprocess.TimeoutExpired:
        log.error("[TIMEOUT-ASSEMBLY] %s 超过 %ss 已中止", task.seq_id, args.task_timeout)
    except FileNotFoundError as exc:
        log.error("[NOT-FOUND] %s: %s\n  若在容器外运行，请确认 %s 已安装；\n"
                  "  或不指定 --mitofinder_script_path 以使用容器内组装器",
                  task.seq_id, exc, args.assembler)
    except Exception as exc:
        log.exception("[EXCEPTION-ASSEMBLY] %s: %s", task.seq_id, exc)
    return False

# ══════════════════════════════════════════════════════════════════════════════
# 样本扫描
# ══════════════════════════════════════════════════════════════════════════════

def collect_assembly_tasks(args: argparse.Namespace) -> List[AssemblyTask]:
    """
    扫描输入目录，返回组装任务列表。
    支持两种方式：
      1. --reads_tsv：TSV 文件（seq_id \\t pe1 \\t pe2 \\t se \\t assembly）
      2. --reads_dir：子目录自动扫描（每个子目录是一个样本）
    """
    tasks: List[AssemblyTask] = []

    # ── 方式 1：TSV 文件 ──────────────────────────────────────────────────────
    if getattr(args, "reads_tsv", None):
        tsv_path = safe_resolve(args.reads_tsv, must_exist=True)
        with tsv_path.open(encoding="utf-8") as f:
            for lineno, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                seq_id = cols[0].strip()

                def _col(idx: int) -> Optional[Path]:
                    if idx < len(cols) and cols[idx].strip() not in ("", "-", "NA"):
                        raw = Path(cols[idx].strip())
                        # [BUG-H] 相对路径基于 TSV 文件所在目录解析，而非运行时 CWD
                        p = raw.resolve() if raw.is_absolute() else (tsv_path.parent / raw).resolve()
                        if not p.exists():
                            log.warning("TSV 行%d: 路径不存在 %s", lineno, p)
                            return None
                        return p
                    return None

                t = AssemblyTask(seq_id=seq_id,
                                 pe1=_col(1), pe2=_col(2), se=_col(3), assembly=_col(4))
                if not any([t.pe1, t.se, t.assembly]):
                    log.warning("TSV 行%d: %s 无有效输入，跳过", lineno, seq_id)
                    continue
                tasks.append(t)
        log.info("从 TSV 读取 %d 个组装任务", len(tasks))
        return tasks

    # ── 方式 2：reads_dir 子目录扫描 ─────────────────────────────────────────
    reads_dir    = safe_resolve(getattr(args, "reads_dir", "."))
    pattern_r1   = getattr(args, "r1_suffix", "_1.clean.fq.gz")
    pattern_r2   = getattr(args, "r2_suffix", "_2.clean.fq.gz")
    pattern_se   = getattr(args, "se_suffix",  ".fastq.gz")
    suffix_subdir = getattr(args, "suffix_subdir", None)

    for sub in sorted(reads_dir.iterdir()):
        if not sub.is_dir():
            continue
        if suffix_subdir and not match_pattern(sub.name, suffix_subdir):
            continue
        sample = sub.name

        r1_files = sorted(sub.glob(f"*{pattern_r1}"))
        r2_files = sorted(sub.glob(f"*{pattern_r2}"))
        if r1_files and r2_files:
            tasks.append(AssemblyTask(seq_id=sample, pe1=r1_files[0].resolve(),
                                      pe2=r2_files[0].resolve()))
            continue

        se_files = [f for f in sorted(sub.glob(f"*{pattern_se}"))
                    if not f.name.endswith(pattern_r1)
                    and not f.name.endswith(pattern_r2)]
        if se_files:
            tasks.append(AssemblyTask(seq_id=sample, se=se_files[0].resolve()))
            continue

        fa_files = sorted(sub.glob("*.fasta")) + sorted(sub.glob("*.fa"))
        if fa_files:
            tasks.append(AssemblyTask(seq_id=sample, assembly=fa_files[0].resolve()))

    if not tasks:
        log.warning("reads_dir=%s 下未扫描到任何样本", reads_dir)

    log.info("扫描到 %d 个样本", len(tasks))
    return tasks


_SPADES_INDICATORS = {"contigs.fasta", "scaffolds.fasta", "assembly_graph.fastg", "spades.log"}


def _safe_sample_id(raw_id: str) -> str:
    """将 FASTA header 第一列转换为可用于文件名和 MitoFinder -j 的样本名。"""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", raw_id).strip("_") or "sample"


def _unique_split_dir(base_dir: Path) -> Path:
    """
    [ADD-MERGED-FASTA] 拆分合并 fasta 时避免覆盖旧文件。
    若 base_dir 已存在，依次使用 base_dir_run2、base_dir_run3 ...
    """
    if not base_dir.exists():
        return ensure_dir(base_dir)
    idx = 2
    while True:
        candidate = base_dir.with_name(f"{base_dir.name}_run{idx}")
        if not candidate.exists():
            return ensure_dir(candidate)
        idx += 1


def _split_merged_fasta_by_id(
    fasta_path: Path,
    args: argparse.Namespace,
    strip_suffix,
) -> List[AnnotateTask]:
    """
    [ADD-MERGED-FASTA] 支持 rename_summary_all.fasta 这类“多样本合并 fasta”。
    规则：以 header 第一列作为样本名分组，例如 ">B17 topology=linear" 归入 B17。
    为避免同一样本内重复 header，写出时改名为 B17_contig1、B17_contig2 ...
    同名 header 不去重、不覆盖，全部作为该样本的多个 contig 保留。

    [OFFICIAL-MITOFINDER] MitoFinder 官方 -a/--assembly 用法支持一个
    assembly.fasta 中包含多条 contig；这里按样本拆分为 B17.fasta、
    B2.fasta 等，正是为了让每个样本的所有 contig 进入同一个注释任务。
    """
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    grouped: Dict[str, List[Tuple[int, SeqRecord]]] = {}
    for rec_idx, rec in enumerate(records, 1):
        seq_id = strip_suffix(_safe_sample_id(rec.id))
        grouped.setdefault(seq_id, []).append((rec_idx, rec))

    split_root = (
        Path(args.output) / "_split_input_fastas"
        if args.output
        else fasta_path.parent / "_split_input_fastas"
    )
    split_dir = _unique_split_dir(split_root / fasta_path.stem)
    manifest_tsv = split_dir / "split_manifest.tsv"

    tasks: List[AnnotateTask] = []
    with manifest_tsv.open("w", newline="", encoding="utf-8") as mf:
        writer = csv.writer(mf, delimiter="\t")
        writer.writerow([
            "sample_name", "new_contig_id", "original_header",
            "original_record_index", "split_fasta",
        ])

        for seq_id, rec_items in sorted(grouped.items()):
            out_fasta = split_dir / f"{seq_id}.fasta"
            out_records: List[SeqRecord] = []
            for idx, (rec_idx, rec) in enumerate(rec_items, 1):
                new_id = f"{seq_id}_contig{idx}"
                desc_tail = (
                    rec.description[len(rec.id):].strip()
                    if rec.description.startswith(rec.id)
                    else rec.description
                )
                description = f"{new_id} {desc_tail}".strip()
                out_records.append(SeqRecord(rec.seq, id=new_id, name=new_id, description=description))
                writer.writerow([seq_id, new_id, rec.description, rec_idx, out_fasta.name])

            SeqIO.write(out_records, str(out_fasta), "fasta")
            tasks.append(AnnotateTask(out_fasta.resolve(), seq_id))

    log.info("合并 fasta 已按 header 第一列拆分: %s → %s  (样本=%d)",
             fasta_path, split_dir, len(tasks))
    log.info("拆分映射表: %s", manifest_tsv)
    # [ADD-MERGED-FASTA] 运行时提示：MitoFinder 官方 -a assembly.fasta
    # 支持一个样本 fasta 内含多条 contig；重复 header 已统一改名，避免覆盖序列。
    log.info(
        "提示: 重复样本名会作为同一样本的多条 contig 保留，并改名为 sample_contigN；"
        "所有原始 header 可在 split_manifest.tsv 中追溯。"
    )
    log.info(
        "MitoFinder 官方 -a/--assembly 支持一个 assembly.fasta 内含一条或多条 contig；"
        "本脚本会将每个样本的 contig 放入同一个样本 fasta 后再注释。"
    )
    return tasks


def _collect_single_fasta_tasks(
    input_path: Path,
    args: argparse.Namespace,
    strip_suffix,
) -> List[AnnotateTask]:
    """
    [ADD-MERGED-FASTA] 单个 fasta 输入的最小兼容层。
    普通单样本 fasta 仍作为 1 个任务；检测到多样本合并 fasta 时才自动拆分。
    """
    records = list(SeqIO.parse(str(input_path), "fasta"))
    ids = [_safe_sample_id(rec.id) for rec in records if rec.id]
    id_counts = {seq_id: ids.count(seq_id) for seq_id in set(ids)}

    # [ADD-MERGED-FASTA] 自动识别当前文件格式：多个样本 ID，且至少一个 ID 出现多次。
    # 这匹配 rename_summary_all.fasta：同一样本可能有多条 contig 记录。
    auto_split = len(id_counts) > 1 and any(n > 1 for n in id_counts.values())
    if getattr(args, "split_input_fasta_by_id", False) or auto_split:
        return _split_merged_fasta_by_id(input_path, args, strip_suffix)

    # [ADD-SINGLE-FASTA] 兼容 annotation2summary -i 直接传入单个 .fasta 文件。
    # 原有目录扫描逻辑保持不变；普通文件输入只在这里生成 1 个注释任务。
    return [AnnotateTask(input_path.resolve(), strip_suffix(input_path.stem))]


def collect_annotate_tasks(args: argparse.Namespace) -> List[AnnotateTask]:
    """扫描 contig fasta，返回注释任务列表。"""
    input_path   = args._input_dir_abs
    suffix_fasta = args.suffix_fasta or ""
    which_fasta  = args.which_fasta2_mitofinder
    tasks: List[AnnotateTask] = []

    def _strip(stem: str) -> str:
        return stem[: -len(suffix_fasta)] if suffix_fasta and stem.endswith(suffix_fasta) else stem

    # [ADD-SINGLE-FASTA][ADD-MERGED-FASTA] 文件输入在这里短路处理；
    # 原有 flat/subdir 目录扫描代码从下面开始，保持不变。
    if input_path.is_file():
        tasks = _collect_single_fasta_tasks(input_path, args, _strip)
        log.info("共找到 %d 个注释任务", len(tasks))
        return tasks

    input_dir = input_path

    if args.fasta_position == "flat":
        for f in sorted(input_dir.glob("*.fasta")):
            if suffix_fasta and suffix_fasta not in f.stem:
                continue
            tasks.append(AnnotateTask(f.resolve(), _strip(f.stem)))

    elif args.fasta_position == "subdir":
        for sub_dir in sorted(input_dir.iterdir()):
            if not sub_dir.is_dir():
                continue
            if not match_pattern(sub_dir.name, args.suffix_subdir):
                continue
            if args.for_spades_dir or which_fasta:
                candidates = [which_fasta] if which_fasta else ["contigs.fasta", "scaffolds.fasta"]
                for cand in candidates:
                    target = sub_dir / cand
                    if target.exists():
                        tasks.append(AnnotateTask(target.resolve(), sub_dir.name))
                        break
                else:
                    log.warning("未找到目标 fasta: %s", sub_dir)
            else:
                for f in sorted(sub_dir.glob("*.fasta")):
                    tasks.append(AnnotateTask(f.resolve(), _strip(f.stem)))

    if not tasks:
        log.warning("未找到符合条件的注释任务，自动诊断 ...")
        ana = _analyze_directory(input_dir)
        log.info("目录结构: %s  推荐参数: %s", ana.structure, ana.recommended)

    log.info("共找到 %d 个注释任务", len(tasks))
    return tasks

# ══════════════════════════════════════════════════════════════════════════════
# 目录智能诊断
# ══════════════════════════════════════════════════════════════════════════════

def _analyze_directory(input_dir: Path) -> DirAnalysis:
    result = DirAnalysis()
    try:
        if not input_dir.exists():
            log.error("目录不存在: %s", input_dir)          # [FIX-5]
            return result
        entries = list(input_dir.iterdir())
    except (PermissionError, FileNotFoundError) as exc:      # [FIX-5]
        log.error("无法读取目录 (%s): %s", type(exc).__name__, input_dir)
        return result

    toplevel_fastas = [f.name for f in sorted(entries) if f.suffix == ".fasta" and f.is_file()]
    subdirs = [d.name for d in sorted(entries) if d.is_dir()]

    if toplevel_fastas:
        result.structure = "flat"
        result.toplevel_fastas = toplevel_fastas[:5]
        result.recommended["fasta_position"] = "flat"
        suffixes = {Path(f).stem.rsplit("_", 1)[-1] for f in toplevel_fastas if "_" in Path(f).stem}
        if len(suffixes) == 1:
            result.recommended["suffix_fasta"] = f"_{suffixes.pop()}"

    if subdirs:
        result.subdirs = subdirs[:5]
        if len(subdirs) > 1:
            parts_list = [d.split("_") for d in subdirs]
            if all(len(p) > 1 for p in parts_list):
                for n in range(1, len(parts_list[0])):
                    candidate = "_".join(parts_list[0][-n:])
                    if all(d.endswith(candidate) for d in subdirs):
                        result.recommended["suffix_subdir"] = f"*{candidate}*"
                        break

        spades_dirs = [s for s in subdirs
                       if any((input_dir / s / ind).exists() for ind in _SPADES_INDICATORS)]
        if spades_dirs:
            result.spades_dirs = spades_dirs
            result.structure = "subdir"
            result.recommended.update({"fasta_position": "subdir", "for_spades_dir": True})
            for cand in ("contigs.fasta", "scaffolds.fasta"):
                if all((input_dir / d / cand).exists() for d in spades_dirs[:3]):
                    result.recommended["which_fasta2_mitofinder"] = cand
                    break
        elif not toplevel_fastas:
            result.structure = "likely_subdir"
            result.recommended["fasta_position"] = "subdir"

    return result

# ══════════════════════════════════════════════════════════════════════════════
# ══  三大主流程  ══════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

# ── A. only_assembly ──────────────────────────────────────────────────────────

def run_only_assembly(args: argparse.Namespace) -> None:
    """
    仅进行批量组装（megahit / metaspades / idba）。
    不调用 MitoFinder，输出原始组装 contig fasta。
    """
    _section("only_assembly：批量独立组装")

    output_dir    = ensure_dir(_resolve_output(args, "assembly_only"))
    setup_logging(output_dir / "assembly_only.log", verbose=args.verbose)
    script_file   = output_dir / "assembly_commands.sh"
    finished_file = output_dir / "assembly_finished.txt"
    finished_set  = load_finished(finished_file)

    tasks = collect_assembly_tasks(args)
    if not tasks:
        log.warning("无组装任务，退出。")
        return

    log.info("并发配置: max_workers=%d × threads/task=%d → 总线程 %d",
             args.max_workers, args.threads, args.max_workers * args.threads)
    if args.max_workers > 1:
        log.info("⚠  注意：并行组装会产生高强度磁盘 I/O，"
                 "请确保 output 目录位于 SSD 或高性能并行文件系统上。")

    ok = fail = 0
    with ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        futures = {
            pool.submit(_run_assembly_only_task, task, output_dir,
                        args, script_file, finished_file, finished_set): task.seq_id
            for task in tasks
        }
        for fut in as_completed(futures):
            exc = fut.exception()
            if exc:
                log.error("任务 %s 异常: %s", futures[fut], exc)
                fail += 1
            elif fut.result():
                ok += 1
            else:
                fail += 1

    log.info("组装完成: 成功=%d  失败=%d", ok, fail)
    _section("only_assembly 完成！")


# ── B. assembly2summary ───────────────────────────────────────────────────────

def run_assembly2summary(args: argparse.Namespace) -> None:
    """
    reads → MitoFinder（组装+注释）→ 提取汇总  全流程。
    MitoFinder 内部调用组装器，无需单独组装步骤。
    """
    _section("assembly2summary：reads → 注释 → 汇总")

    output_dir    = ensure_dir(_resolve_output(args, "mitofinder_assembly"))
    setup_logging(output_dir / "assembly2summary.log", verbose=args.verbose)
    script_file   = output_dir / "mf_assembly_scripts.sh"
    finished_file = output_dir / "mf_assembly_finished.txt"
    finished_set  = load_finished(finished_file)

    tasks = collect_assembly_tasks(args)
    if not tasks:
        log.warning("无组装任务，退出。")
        return
    # [ADD-RUNTIME-HINT] 非 dry-run 时先确认 MitoFinder 可执行环境；
    # 避免批量样本全部启动后才逐个报 mitofinder not found。
    if not args.dry_run and not _mitofinder_runtime_ready(args):
        _log_mitofinder_runtime_hint()
        return

    genbank_ref = args._genbank_ref_abs

    def _iter():
        for task in tasks:
            mf_args    = _build_mf_args_assembly(task, genbank_ref, args)
            bind_paths = [p for p in [task.pe1, task.pe2, task.se, task.assembly, genbank_ref] if p]
            yield task.seq_id, mf_args, bind_paths

    _run_pool(_iter(), output_dir, args, script_file, finished_file, finished_set)

    _section("assembly2summary：开始提取汇总")
    run_summary(args, annotation_dir=output_dir)
    _section("assembly2summary 完成！")


# ── C. annotation2summary ─────────────────────────────────────────────────────

def run_annotation2summary(args: argparse.Namespace) -> None:
    """
    已有 contig fasta → MitoFinder 注释 → 提取汇总。
    """
    _section("annotation2summary：注释 → 汇总")

    output_dir    = ensure_dir(_resolve_output(args, "mitofinder_annotate"))
    setup_logging(output_dir / "annotation2summary.log", verbose=args.verbose)
    script_file   = output_dir / "mf_annotate_scripts.sh"
    finished_file = output_dir / "mf_annotate_finished.txt"
    finished_set  = load_finished(finished_file)

    tasks = collect_annotate_tasks(args)
    if not tasks:
        log.warning("无注释任务，退出。")
        return
    # [ADD-RUNTIME-HINT] 非 dry-run 时先确认 MitoFinder 可执行环境；
    # 避免批量样本全部启动后才逐个报 mitofinder not found。
    if not args.dry_run and not _mitofinder_runtime_ready(args):
        _log_mitofinder_runtime_hint()
        return

    genbank_ref = args._genbank_ref_abs

    def _iter():
        for task in tasks:
            mf_args    = _build_mf_args_annotate(task.seq_id, task.fasta_path, genbank_ref, args)
            bind_paths = [task.fasta_path, genbank_ref]
            yield task.seq_id, mf_args, bind_paths

    _run_pool(_iter(), output_dir, args, script_file, finished_file, finished_set)

    _section("annotation2summary：开始提取汇总")
    run_summary(args, annotation_dir=output_dir)
    _section("annotation2summary 完成！")

# ══════════════════════════════════════════════════════════════════════════════
# 提取 / 汇总（被三大流程共用）
# ══════════════════════════════════════════════════════════════════════════════

def run_summary(args: argparse.Namespace, annotation_dir: Optional[Path] = None) -> None:
    """对注释结果目录执行全套提取和汇总操作。"""
    if annotation_dir is None:
        annotation_dir = safe_resolve(args.extract_input_dir)

    summary_dir = ensure_dir(
        Path(args.output) / "mitofinder_summary" if args.output
        else safe_resolve(args.extract_output_dir)
    )
    setup_logging(summary_dir / "summary.log", verbose=args.verbose)

    sfx     = args.suffix2_mitofinder_results
    rm_sfx  = args.if_write_basename2_summary_by_remove_suffix
    sub_pat = args.suffix_subdir

    # 1. 基因长度汇总表
    summary_basic_information(
        annotation_dir,
        summary_dir / "summary_mitofinder_results.tsv",
        sub_pat,
    )
    # 2. 按基因收集序列
    gather_genes(annotation_dir, summary_dir / "gather_genes", args.genes, sub_pat)

    # 3. 合并 mtDNA contigs
    gather_mt_contigs(annotation_dir, summary_dir / "mt_contigs",
                      args.extension2_contigs, sub_pat)

    # 4. 收集 contig fasta / gb
    collect_assembly_results(
        annotation_dir,
        summary_dir / args.output2_mitofinder_assembly_result,
        sfx, rm_sfx,
    )

    # 5. 汇总 .infos（可选）
    if args.summarize_assembly_infos:
        summarize_assembly_infos(
            annotation_dir,
            summary_dir / args.output2_mitofinder_assembly_infos_result,
            sfx, rm_sfx,
        )

    # 6. 汇总 tRNA（可选）
    if args.summarize_trna:
        summarize_trna_results(
            annotation_dir,
            ensure_dir(summary_dir / args.output2_mitofinder_trna_result),
            sfx, rm_sfx,
        )

    log.info("汇总输出目录: %s", summary_dir)

# ══════════════════════════════════════════════════════════════════════════════
# Snakemake 工作流生成器
# ══════════════════════════════════════════════════════════════════════════════

_SNAKEFILE_TEMPLATE = '''\
# =============================================================================
# MitoFinder Snakemake Workflow  (auto-generated by mitofinder_pipeline.py v{version})
# 用法（本地）：
#   snakemake -s Snakefile --configfile config.yaml --cores {total_cores}
# 用法（SLURM）：
#   snakemake -s Snakefile --configfile config.yaml --jobs {max_workers} \\
#             --cluster "sbatch --ntasks={threads} --mem={memory}G --time=12:00:00"
# =============================================================================

configfile: "config.yaml"

SIF          = config["sif_path"]
GENBANK_REF  = config["genbank_ref"]
OUTPUT_DIR   = config["output_dir"]
SAMPLES      = config["samples"]
THREADS      = config["threads"]
MEMORY       = config["memory"]
GENETIC_CODE = config["genetic_code"]
TOOL_PARAM   = config["tool_param"]
ASSEMBLER    = config.get("assembler", "megahit")
SNAKE_MODE   = config.get("snake_mode", "annotate")   # assembly | annotate


rule all:
    input:
        expand(OUTPUT_DIR + "/{{sample}}/{{sample}}_Final_Results", sample=SAMPLES)


rule mitofinder:
    input:
        reads  = lambda wc: _get_input_files(wc.sample),
        refseq = GENBANK_REF,
    output:
        directory(OUTPUT_DIR + "/{{sample}}/{{sample}}_Final_Results")
    threads: THREADS
    resources:
        mem_gb = MEMORY
    params:
        outdir = OUTPUT_DIR,
        extra  = config.get("extra_flags", ""),
    log:
        OUTPUT_DIR + "/logs/{{sample}}.log"
    run:
        import os, subprocess
        info       = SAMPLES[wildcards.sample]
        reads_flag = _reads_flags(info)
        ref_dir    = os.path.dirname(os.path.abspath(input.refseq))
        out_dir    = os.path.abspath(params.outdir)

        # [BUG-B] --bind 路径去重：ref/reads 同目录时避免 FATAL: Duplicate bind path
        reads_files = input.reads if input.reads else []
        reads_dir   = os.path.dirname(os.path.abspath(reads_files[0])) if reads_files else ref_dir
        bind_paths  = list(dict.fromkeys([out_dir, ref_dir, reads_dir]))  # 保序去重
        bind_str    = " ".join(f"--bind {{b}}:{{b}}" for b in bind_paths if b)

        cmd = (
            f"singularity run "
            f"{{bind_str}} "
            f"--pwd {{out_dir}} "
            f"{{SIF}} "
            f"{{reads_flag}} "
            f"-j {{wildcards.sample}} "
            f"-r {{input.refseq}} "
            f"-t {{TOOL_PARAM}} "
            f"-o {{GENETIC_CODE}} "
            f"-p {{threads}} "
            f"-m {{resources.mem_gb}} "
            f"--adjust-direction "
            f"{{params.extra}} "
            f"> {{log}} 2>&1"
        )
        subprocess.run(cmd, shell=True, check=True, cwd=out_dir)


def _get_input_files(sample):
    info = SAMPLES[sample]
    return [v for k in ("pe1", "pe2", "se", "assembly")
            for v in ([info[k]] if info.get(k, "") not in ("", "NA", "-") else [])]


def _reads_flags(info):
    """
    [BUG-B2] 仅在 reads 模式下追加组装器 Flag（--megahit 等）。
    已有 assembly 时只传 -a，绝不追加组装器 Flag，避免 MitoFinder 参数冲突。
    """
    if info.get("pe1") and info.get("pe2"):
        return f"-1 {{info[\'pe1\']}} -2 {{info[\'pe2\']}} --{{ASSEMBLER}}"
    elif info.get("se"):
        return f"-s {{info[\'se\']}} --{{ASSEMBLER}}"
    elif info.get("assembly"):
        return f"-a {{info[\'assembly\']}}"   # 注释模式，不附加组装器 Flag
    return ""
'''

_CONFIG_TEMPLATE = """\
# MitoFinder Snakemake config (auto-generated by mitofinder_pipeline.py v{version})
sif_path:     "{sif_path}"
genbank_ref:  "{genbank_ref}"
output_dir:   "{output_dir}"
threads:      {threads}
memory:       {memory}
genetic_code: {genetic_code}
tool_param:   "{tool_param}"
assembler:    "{assembler}"
snake_mode:   "{snake_mode}"
extra_flags:  ""

# samples: seq_id -> {{pe1, pe2, se, assembly}}  (空字段填 "")
samples:
{samples_yaml}
"""


def run_snakemake(args: argparse.Namespace) -> None:
    """生成 Snakefile + config.yaml + run_snakemake.sh。"""
    _section("snakemake：生成工作流文件")

    out = ensure_dir(
        Path(args.output) if args.output else Path("mitofinder_snakemake")
    )
    setup_logging(out / "snakemake_gen.log", verbose=args.verbose)

    snake_mode = getattr(args, "snake_mode", "annotate")
    sif_path   = str(safe_resolve(args.singularity_sif)) if args.singularity_sif else _DEFAULT_SIF
    genbank_ref = str(args._genbank_ref_abs) if args._genbank_ref_abs else ""
    output_dir  = str(out / "results")
    assembler   = getattr(args, "assembler", "megahit")
    total_cores = args.max_workers * args.threads

    # 收集样本
    if snake_mode == "assembly":
        raw_tasks = collect_assembly_tasks(args)
        yaml_lines = []
        for t in raw_tasks:
            yaml_lines += [
                f"  {t.seq_id}:",
                f'    pe1:      "{t.pe1 or ""}"',
                f'    pe2:      "{t.pe2 or ""}"',
                f'    se:       "{t.se or ""}"',
                f'    assembly: "{t.assembly or ""}"',
            ]
    else:
        raw_tasks = collect_annotate_tasks(args)
        yaml_lines = []
        for t in raw_tasks:
            yaml_lines += [
                f"  {t.seq_id}:",
                f'    pe1:      ""',
                f'    pe2:      ""',
                f'    se:       ""',
                f'    assembly: "{t.fasta_path}"',
            ]

    config_txt = _CONFIG_TEMPLATE.format(
        version=__version__,
        sif_path=sif_path, genbank_ref=genbank_ref,
        output_dir=output_dir, threads=args.threads,
        memory=args.memory, genetic_code=args.genetic_code,
        tool_param=args.tool_param, assembler=assembler,
        snake_mode=snake_mode,
        samples_yaml="\n".join(yaml_lines),
    )
    snake_txt = _SNAKEFILE_TEMPLATE.format(
        version=__version__,
        total_cores=total_cores, threads=args.threads,
        memory=args.memory, max_workers=args.max_workers,
    )
    run_sh = textwrap.dedent(f"""\
        #!/bin/bash
        # ── 本地多核运行 ─────────────────────────────────────────────────────────
        # 注意：Snakefile 已通过 subprocess 手动调用 singularity run，
        #       无需 --use-singularity（重复加此参数在某些集群会引发警告或报错）
        snakemake -s {out}/Snakefile \\
          --configfile {out}/config.yaml \\
          --cores {total_cores} \\
          --printshellcmds --reason

        # ── HPC SLURM 投递（取消注释后使用）─────────────────────────────────────
        # snakemake -s {out}/Snakefile \\
        #   --configfile {out}/config.yaml \\
        #   --jobs {args.max_workers} \\
        #   --cluster "sbatch --ntasks={args.threads} --mem={args.memory}G --time=12:00:00"
    """)

    (out / "Snakefile").write_text(snake_txt, encoding="utf-8")
    (out / "config.yaml").write_text(config_txt, encoding="utf-8")
    (out / "run_snakemake.sh").write_text(run_sh, encoding="utf-8")
    (out / "run_snakemake.sh").chmod(0o755)

    log.info("Snakefile      : %s", out / "Snakefile")
    log.info("config.yaml    : %s", out / "config.yaml")
    log.info("run_snakemake.sh: %s", out / "run_snakemake.sh")
    log.info("样本数: %d", len(raw_tasks))
    log.info("执行: bash %s", out / "run_snakemake.sh")
    _section("snakemake 生成完成！")

# ══════════════════════════════════════════════════════════════════════════════
# 提取子功能实现（汇总用）
# ══════════════════════════════════════════════════════════════════════════════

def _parse_nt_fasta(fasta_path: Path) -> Dict[str, str]:
    """[FIX-3] 重复基因长度逗号拼接。"""
    gene_lengths: Dict[str, str] = {}
    try:
        for rec in SeqIO.parse(str(fasta_path), "fasta"):
            if "@" in rec.description:
                gene = rec.description.split("@", 1)[1].strip()
                l = str(len(rec.seq))
                gene_lengths[gene] = f"{gene_lengths[gene]},{l}" if gene in gene_lengths else l
    except Exception as exc:
        log.warning("解析 NT fasta 失败 %s: %s", fasta_path, exc)
    return gene_lengths


def summary_basic_information(
    input_dir: Path, output_tsv: Path,
    suffix_subdir: Optional[str] = None,
) -> None:
    ensure_dir(output_tsv.parent)
    sample_data: Dict[str, Dict[str, str]] = {}
    all_genes: Set[str] = set()

    for _, final_dir, sample_name in iter_final_results_dirs(input_dir, suffix_subdir):
        for nt_fasta in sorted(final_dir.glob("*_final_genes_NT.fasta")):
            genes = _parse_nt_fasta(nt_fasta)
            if genes:
                sample_data[sample_name] = genes
                all_genes.update(genes)

    sorted_genes = sorted(all_genes)
    with output_tsv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Sample_Name"] + sorted_genes)
        for sample, genes in sorted(sample_data.items()):
            writer.writerow([sample] + [genes.get(g, "N") for g in sorted_genes])

    log.info("基因信息汇总 → %s  (样本=%d, 基因=%d)",
             output_tsv, len(sample_data), len(sorted_genes))


def gather_genes(
    input_dir: Path, output_dir: Path,
    genes: Optional[List[str]] = None,
    suffix_subdir: Optional[str] = None,
) -> None:
    ensure_dir(output_dir)
    if genes is None:
        gene_set: Set[str] = set()
        for _, final_dir, _ in iter_final_results_dirs(input_dir, suffix_subdir):
            for nt_fasta in final_dir.glob("*_final_genes_NT.fasta"):
                gene_set.update(_parse_nt_fasta(nt_fasta).keys())
        genes = sorted(gene_set)
        log.info("自动识别 %d 个基因: %s", len(genes), genes)

    if not genes:
        log.warning("未发现任何基因，gather_genes 跳过。")
        return

    gene_records: Dict[str, List] = {g: [] for g in genes}
    for _, final_dir, _ in iter_final_results_dirs(input_dir, suffix_subdir):
        for nt_fasta in sorted(final_dir.glob("*_final_genes_NT.fasta")):
            try:
                for rec in SeqIO.parse(str(nt_fasta), "fasta"):
                    if "@" in rec.description:
                        gname = rec.description.split("@", 1)[1].strip()
                        if gname in gene_records:
                            existing = gene_records[gname]
                            # [BUG-J] 同名基因多拷贝时追加 _copyN 后缀保证 ID 唯一，
                            # 避免 MAFFT 等 MSA 工具因重复 Header 而报错中断
                            dup_count = sum(1 for r in existing if r.id == rec.id)
                            if dup_count > 0:
                                rec = rec.__class__(
                                    rec.seq,
                                    id=f"{rec.id}_copy{dup_count + 1}",
                                    description=rec.description,
                                )
                            existing.append(rec)
            except Exception as exc:
                log.warning("读取 %s 失败: %s", nt_fasta, exc)

    for gene, recs in gene_records.items():
        try:
            SeqIO.write(recs, str(output_dir / f"{gene}.fasta"), "fasta")
            log.info("  %s: %d 条", gene, len(recs))
        except Exception as exc:
            log.warning("写出 %s 失败: %s", gene, exc)


def gather_mt_contigs(                                       # [FIX-2]
    input_dir: Path, output_dir: Path,
    extensions: str = ".fasta",
    suffix_subdir: Optional[str] = None,
) -> None:
    ensure_dir(output_dir)
    exts = [e.strip() for e in extensions.split(",") if e.strip()]
    ext_pat = "|".join(re.escape(e) for e in exts)
    contig_re = re.compile(rf".*_mtDNA_contig(?:_\d+)?({ext_pat})$")   # (?:_\d+)? 可选

    for _, final_dir, sample_name in iter_final_results_dirs(input_dir, suffix_subdir):
        recs: List[SeqRecord] = []
        for f in sorted(final_dir.iterdir()):
            if f.is_file() and contig_re.match(f.name):
                try:
                    recs.extend(SeqIO.parse(str(f), "fasta"))
                except Exception as exc:
                    log.warning("读取 %s 失败: %s", f.name, exc)
        if not recs:
            log.warning("样本 %s: 未找到 mtDNA contig", sample_name)
            continue
        try:
            SeqIO.write(recs, str(output_dir / f"{sample_name}_all_mt.fasta"), "fasta")
            with (output_dir / f"{sample_name}_mt_contigsID_list.txt").open("w") as fh:
                for rec in recs:
                    fh.write(re.sub(r"\(.*?\)", "", rec.description).strip() + "\n")
            log.info("  %s: %d 条 contig", sample_name, len(recs))
        except Exception as exc:
            log.warning("写出 %s 失败: %s", sample_name, exc)


_FASTA_PAT = re.compile(r".*_mtDNA_contig(?:_\d+)?\.fasta$")
_GB_PAT    = re.compile(r".*_mtDNA_contig(?:_\d+)?\.gb$")


def collect_assembly_results(                               # [FIX-4]
    input_dir: Path, output_dir: Path,
    suffix: str = "", remove_suffix: str = "no",
) -> Tuple[int, int]:
    fasta_dir = ensure_dir(output_dir / "fasta_files")
    gb_dir    = ensure_dir(output_dir / "gb_files")
    n_fasta = n_gb = 0

    def _name(raw: str) -> str:
        return raw[: -len(suffix)] if (remove_suffix.lower() == "yes" and suffix and raw.endswith(suffix)) else raw

    for _, final_dir, raw_name in iter_final_results_dirs(input_dir):
        name = _name(raw_name)
        for f in sorted(final_dir.iterdir()):
            if not f.is_file():
                continue
            if _FASTA_PAT.match(f.name):
                tdir, n_fasta = fasta_dir, n_fasta + 1
            elif _GB_PAT.match(f.name):
                tdir, n_gb = gb_dir, n_gb + 1
            else:
                continue
            # [BUG-D] 严格边界匹配：startswith(name + "_") 避免 S1 误匹配 S10_xxx.fasta
            new_name = f.name if f.name.startswith(f"{name}_") else f"{name}_{f.name}"
            try:
                shutil.copy2(f, tdir / new_name)
            except Exception as exc:
                log.warning("复制 %s 失败: %s", f.name, exc)

    (output_dir / "assembly_results_summary.txt").write_text(
        f"FASTA={n_fasta}  GB={n_gb}\nfasta_dir: {fasta_dir}\ngb_dir: {gb_dir}\n",
        encoding="utf-8",
    )
    log.info("收集完成: FASTA=%d, GB=%d", n_fasta, n_gb)
    return n_fasta, n_gb


_INFOS_FIELDS = {
    "statistics":      re.compile(r"Statistics for\s*[:：]?\s*(.+)", re.I),
    "initial_contig":  re.compile(r"Initial contig name\s*[:：]\s*(.+)"),
    "length":          re.compile(r"Length\s*[:：]\s*([\d\.]+)"),
    "gc_content":      re.compile(r"GC content\s*[:：]\s*([\d\.]+)%"),
    "circularization": re.compile(r"Circularization\s*[:：]\s*(.+)"),
}


def summarize_assembly_infos(
    input_dir: Path, output_tsv: Path,
    suffix: str = "", remove_suffix: str = "no",
) -> None:
    ensure_dir(output_tsv.parent)

    def _name(raw: str) -> str:
        return raw[: -len(suffix)] if (remove_suffix.lower() == "yes" and suffix and raw.endswith(suffix)) else raw

    rows: List[list] = []
    for _, final_dir, raw_name in iter_final_results_dirs(input_dir):
        name = _name(raw_name)
        infos = sorted(final_dir.glob("*.infos"))
        if not infos:
            rows.append([name, "", "", "", "", ""])
            continue
        for inf in infos:
            try:
                text = inf.read_text(encoding="utf-8", errors="replace")
                row  = {k: (m.group(1).strip() if (m := pat.search(text)) else "")
                        for k, pat in _INFOS_FIELDS.items()}
                stats = row["statistics"]
                if len(infos) == 1 and stats.lower() == "final sequence":
                    stats = "contig1"
                rows.append([name, stats.replace(":", "").replace(" ", "_"),
                             row["initial_contig"], row["length"],
                             row["gc_content"], row["circularization"]])
            except Exception as exc:
                log.warning("解析 %s 失败: %s", inf, exc)
                rows.append([name, "", "", "", "", ""])

    with output_tsv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Result_Dir", "Statistics", "Initial_contig_name",
                         "Length", "GC_content", "Circularization"])
        writer.writerows(rows)
    log.info("Assembly infos → %s (%d 行)", output_tsv, len(rows))


_TRNA_KEYWORDS  = ("_trnascan", "_arwen", "_mitfi")
_TRNA_HEADER_RE = re.compile(r"(Sequence|Name)\s+t(RNA|RNA #)\s+Bounds?\s+tRNA\s+Anti", re.I)


def summarize_trna_results(
    input_dir: Path, output_dir: Path,
    suffix: str = "", remove_suffix: str = "no",
) -> None:
    ensure_dir(output_dir)

    def _name(raw: str) -> str:
        return raw[: -len(suffix)] if (remove_suffix.lower() == "yes" and suffix and raw.endswith(suffix)) else raw

    for _, final_dir, raw_name in iter_final_results_dirs(input_dir):
        name = _name(raw_name)
        for trna_dir in sorted(final_dir.iterdir()):
            if not (trna_dir.is_dir() and any(trna_dir.name.endswith(k) for k in _TRNA_KEYWORDS)):
                continue
            headers: List[str] = []
            data:    List[str] = []
            try:
                for tf in sorted(trna_dir.iterdir()):
                    if not tf.is_file() or ".log" in tf.name.lower():
                        continue
                    lines = tf.read_text(encoding="utf-8", errors="replace").splitlines()
                    in_data = False
                    for i, line in enumerate(lines):
                        if _TRNA_HEADER_RE.search(line):
                            in_data = True
                            if not headers:
                                headers.append(line)
                                if i + 1 < len(lines) and lines[i + 1].startswith("---"):
                                    headers.append(lines[i + 1])
                        elif in_data and line.strip() and not line.startswith(("Name", "Sequence", "---")):
                            if "\t" in line or "  " in line:
                                data.append(line)
            except Exception as exc:
                log.warning("读取 tRNA %s 失败: %s", trna_dir, exc)
                continue
            if headers and data:
                out_csv = output_dir / f"{name}_{trna_dir.name}_summary.csv"
                with out_csv.open("w", encoding="utf-8") as fh:
                    fh.writelines(h + "\n" for h in headers)
                    fh.writelines(d + "\n" for d in data)
                log.info("  tRNA → %s (%d 行)", out_csv.name, len(data))

# ══════════════════════════════════════════════════════════════════════════════
# 辅助：输出目录解析
# ══════════════════════════════════════════════════════════════════════════════

def _resolve_output(args: argparse.Namespace, default_subdir: str) -> Path:
    """根据 --output 参数决定输出目录。"""
    if args.output:
        return Path(args.output) / default_subdir
    return safe_resolve(getattr(args, f"{default_subdir}_dir", default_subdir))

# ══════════════════════════════════════════════════════════════════════════════
# 命令行参数
# ══════════════════════════════════════════════════════════════════════════════

def _add_common_args(p: argparse.ArgumentParser) -> None:
    """向子命令解析器添加所有模式共用的参数。"""
    # 全局
    p.add_argument("--output", "-o", default=None, metavar="DIR",
                   help="总输出根目录（绝对或相对均可）")
    p.add_argument("--verbose",      action="store_true", help="显示 DEBUG 日志")
    p.add_argument("--dry_run",      action="store_true", help="打印命令但不执行")

    # MitoFinder 核心
    mf = p.add_argument_group("MitoFinder 核心")
    mf.add_argument("--genbank_ref", "-r", metavar="FILE",
                    help="参考线粒体 GenBank 文件（.gb）")
    mf.add_argument("--genetic_code", type=int, default=5, metavar="N",
                    help="NCBI 遗传密码编号（5=无脊椎动物线粒体）")
    mf.add_argument("--threads",  "-p", type=int, default=4,   metavar="N",
                    help="每任务 CPU 线程数")
    mf.add_argument("--memory",   "-m", type=int, default=50,  metavar="GB",
                    help="单任务内存上限（GB）")
    mf.add_argument("--tool_param", "-t", default="mitfi",
                    choices=["mitfi", "trnascan", "arwen"],
                    help="tRNA 注释工具")
    mf.add_argument("--singularity_sif", default=None, metavar="FILE",
                    help=f"Singularity SIF 路径（默认 {_DEFAULT_SIF}）")
    mf.add_argument("--mitofinder_script_path", default="", metavar="PATH",
                    help="MitoFinder 可执行路径；留空则走 Singularity")

    # MitoFinder 高级
    adv = p.add_argument_group("MitoFinder 高级参数")
    adv.add_argument("--adjust_direction", action="store_true", default=True)
    adv.add_argument("--override",         action="store_true")
    adv.add_argument("--ignore_non_std",   action="store_true")
    adv.add_argument("--new_genes",        action="store_true")
    adv.add_argument("--allow_intron",     action="store_true")
    adv.add_argument("--out_gb",           action="store_true")
    adv.add_argument("--cds_merge",        action="store_true")
    adv.add_argument("--blast_eval",        type=float, default=None, metavar="FLOAT")
    adv.add_argument("--nwalk",             type=int,   default=None, metavar="N")
    adv.add_argument("--min_contig_size",   type=int,   default=None, metavar="BP")
    adv.add_argument("--max_contig_size",   type=int,   default=None, metavar="BP")
    adv.add_argument("--circular_offset",   type=int,   default=None, metavar="BP",
                     help="环化偏移量（默认 200）")
    adv.add_argument("--circular_size",     type=int,   default=None, metavar="N",
                     help="环化大小（默认 45）")
    adv.add_argument("--blast_identity_nucl", type=float, default=None, metavar="PCT")
    adv.add_argument("--max_contig",        type=int,   default=None, metavar="N")
    adv.add_argument("--intron_size",       type=int,   default=None, metavar="BP")
    adv.add_argument("--rename_contig",     choices=["yes", "no"], default=None)
    adv.add_argument("--config_file",       default=None, metavar="FILE")

    # 并发
    cc = p.add_argument_group("并发控制")
    cc.add_argument("--max_workers", type=int, default=4,  metavar="N",
                    help="并行任务数（建议 ≈ 总CPU / threads）")
    cc.add_argument("--task_timeout", type=int, default=None, metavar="SEC",
                    help="单任务超时秒数")

    # reads 输入
    ri = p.add_argument_group("reads 输入")
    ri.add_argument("--reads_dir",  default=".", metavar="DIR",
                    help="含样本子目录的 reads 根目录")
    ri.add_argument("--reads_tsv",  default=None, metavar="FILE",
                    help="TSV 文件（seq_id\\tpe1\\tpe2\\tse\\tassembly）")
    ri.add_argument("--r1_suffix",  default="_1.clean.fq.gz", metavar="STR",
                    help="forward reads 文件名后缀")
    ri.add_argument("--r2_suffix",  default="_2.clean.fq.gz", metavar="STR",
                    help="reverse reads 文件名后缀")
    ri.add_argument("--se_suffix",  default=".fastq.gz",      metavar="STR",
                    help="single-end reads 后缀")
    ri.add_argument("--assembler",  default="megahit",
                    choices=["megahit", "metaspades", "idba"],
                    help="组装器")
    ri.add_argument("--shortest_contig", type=int, default=None, metavar="BP")
    ri.add_argument("--suffix_subdir",   default=None,
                    help="fnmatch 模式过滤样本子目录")

    # contig 注释输入
    an = p.add_argument_group("contig 注释输入")
    # [ADD-SINGLE-FASTA] 仅更新帮助文本，说明 -i 现在也可直接接收单个 .fasta 文件。
    an.add_argument("--input_dir2_annotate", "-i", default=".", metavar="PATH",
                    help="含 contig fasta 的输入目录，或单个 fasta 文件")
    # [ADD-MERGED-FASTA] 手动强制按 header 第一列拆分单个合并 fasta；
    # 默认也会自动识别 rename_summary_all.fasta 这类含重复样本 ID 的合并文件。
    an.add_argument("--split_input_fasta_by_id", action="store_true",
                    help="单个 fasta 输入时，按 header 第一列拆成多个样本任务")
    an.add_argument("--fasta_position", choices=["flat", "subdir"], default="flat")
    an.add_argument("--suffix_fasta",   default=None)
    an.add_argument("--for_spades_dir", action="store_true")
    an.add_argument("--which_fasta2_mitofinder", default=None, metavar="FILE")

    # 汇总参数
    sm = p.add_argument_group("汇总参数")
    sm.add_argument("--extract_input_dir",  default="mitofinder_annotate", metavar="DIR",
                    help="提取步骤输入目录（仅 extract_only 模式使用）")
    sm.add_argument("--extract_output_dir", default="mitofinder_summary",  metavar="DIR",
                    help="提取步骤输出目录（仅 extract_only 模式使用）")
    sm.add_argument("--genes", "-g", nargs="*", default=None, metavar="GENE",
                    help="指定基因名；不指定则自动发现")
    sm.add_argument("--extension2_contigs", default=".fasta",
                    help="mtDNA contig 扩展名，逗号分隔")
    sm.add_argument("--summarize_assembly_infos", action="store_true")
    sm.add_argument("--summarize_trna",           action="store_true")
    sm.add_argument("--output2_mitofinder_assembly_result",
                    default="mitofinder_assembly_results")
    sm.add_argument("--output2_mitofinder_assembly_infos_result",
                    default="mitofinder_assembly_infos_summary.tsv")
    sm.add_argument("--output2_mitofinder_trna_result",
                    default="mitofinder_trna_summary")
    sm.add_argument("--suffix2_mitofinder_results", default="")
    sm.add_argument("--if_write_basename2_summary_by_remove_suffix",
                    default="no", choices=["yes", "no"])


def build_parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(
        prog="mitofinder_pipeline",
        description="MitoFinder 整合处理工具 v" + __version__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
        示例：
          # 仅组装
          %(prog)s only_assembly --reads_dir /data/reads -o /data/out -r ref.gb

          # reads → 注释 → 汇总
          %(prog)s assembly2summary --reads_dir /data/reads -o /data/out -r ref.gb

          # 已有 contig 目录 → 注释 → 汇总
          %(prog)s annotation2summary -i /data/contigs -o /data/out -r ref.gb

          # 单个 fasta → 注释 → 汇总
          %(prog)s annotation2summary -i /data/contigs/sample.fasta -o /data/out -r ref.gb

          # 生成 Snakemake 工作流
          %(prog)s snakemake --reads_dir /data/reads -o /data/wf -r ref.gb
        """),
    )
    root.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    sub = root.add_subparsers(dest="subcmd", metavar="SUBCOMMAND")
    sub.required = True

    # ── only_assembly ─────────────────────────────────────────────────────────
    p_oa = sub.add_parser(
        "only_assembly",
        help="仅批量组装 reads（megahit / metaspades / idba），不做 MitoFinder 注释",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_common_args(p_oa)

    # ── assembly2summary ──────────────────────────────────────────────────────
    p_a2s = sub.add_parser(
        "assembly2summary",
        help="reads → MitoFinder 组装+注释 → 提取汇总  全流程",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_common_args(p_a2s)

    # ── annotation2summary ────────────────────────────────────────────────────
    p_n2s = sub.add_parser(
        "annotation2summary",
        help="已有 contig fasta → MitoFinder 注释 → 提取汇总",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_common_args(p_n2s)

    # ── snakemake ─────────────────────────────────────────────────────────────
    p_snk = sub.add_parser(
        "snakemake",
        help="生成 Snakemake 工作流文件（支持 HPC 集群投递）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_common_args(p_snk)
    p_snk.add_argument("--snake_mode", choices=["assembly", "annotate"], default="annotate",
                       help="Snakemake 工作流类型")

    # ── analyze ───────────────────────────────────────────────────────────────
    p_ana = sub.add_parser(
        "analyze",
        help="智能分析输入目录结构，给出参数建议",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p_ana.add_argument("input_dir", metavar="DIR", help="要分析的目录")
    p_ana.add_argument("--verbose", action="store_true")

    return root

# ══════════════════════════════════════════════════════════════════════════════
# 主入口
# ══════════════════════════════════════════════════════════════════════════════

BANNER = r"""
╔══════════════════════════════════════════════════════════════════════╗
║         MitoFinder 整合处理工具 v2.6                                ║
╠══════════════════════════════════════════════════════════════════════╣
║  only_assembly      仅批量组装 reads（不做注释）                     ║
║  assembly2summary   reads → 组装+注释 → 汇总  全流程                ║
║  annotation2summary 已有 contig → 注释 → 汇总                       ║
║  snakemake          生成 Snakemake 工作流（支持 HPC）                ║
║  analyze            智能分析输入目录结构                              ║
╚══════════════════════════════════════════════════════════════════════╝"""


def main() -> None:
    print(BANNER)
    parser = build_parser()
    args   = parser.parse_args()

    setup_logging(verbose=getattr(args, "verbose", False))

    # ── analyze 子命令（无需路径解析）─────────────────────────────────────────
    if args.subcmd == "analyze":
        d = safe_resolve(args.input_dir)
        log.info("分析目录: %s", d)
        ana = _analyze_directory(d)
        print(f"\n目录结构: {ana.structure}")
        if ana.toplevel_fastas: print("顶层 fasta 示例:", ana.toplevel_fastas)
        if ana.subdirs:         print("子目录示例:",      ana.subdirs)
        if ana.spades_dirs:     print("SPAdes 目录:",     ana.spades_dirs)
        if ana.recommended:
            print("\n推荐参数:")
            for k, v in ana.recommended.items():
                print(f"  --{k} {v}")
        return

    # ── 统一路径解析 ──────────────────────────────────────────────────────────
    args._input_dir_abs   = safe_resolve(getattr(args, "input_dir2_annotate", "."))
    args._genbank_ref_abs = safe_resolve(args.genbank_ref) if args.genbank_ref else None
    if args.output:
        args.output = str(safe_resolve(args.output))

    # ── 参数校验 ──────────────────────────────────────────────────────────────
    if args.subcmd in ("assembly2summary", "annotation2summary", "snakemake"):
        if not args.genbank_ref:
            parser.error(f"子命令 {args.subcmd} 需要 --genbank_ref 参数。")

    if args.subcmd == "annotation2summary":
        if not args._input_dir_abs.exists():
            parser.error(f"注释输入路径不存在: {args._input_dir_abs}")
        # [ADD-SINGLE-FASTA] 原代码只允许 -i 是目录；这里最小放宽为“目录或单个 .fasta 文件”。
        # 目录输入仍走原来的 collect_annotate_tasks 扫描逻辑，不改变原有批量目录用法。
        if args._input_dir_abs.is_file() and args._input_dir_abs.suffix != ".fasta":
            parser.error(f"注释输入文件不是 .fasta 后缀: {args._input_dir_abs}")
        if not (args._input_dir_abs.is_dir() or args._input_dir_abs.is_file()):
            parser.error(f"注释输入路径不是文件或目录: {args._input_dir_abs}")

    # ── 分发 ──────────────────────────────────────────────────────────────────
    dispatch = {
        "only_assembly":      run_only_assembly,
        "assembly2summary":   run_assembly2summary,
        "annotation2summary": run_annotation2summary,
        "snakemake":          run_snakemake,
    }
    dispatch[args.subcmd](args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[中断] 程序被用户终止。")
    except Exception as exc:
        log.exception("程序异常退出: %s", exc)
        if "--debug" in sys.argv:
            raise
        sys.exit(1)
