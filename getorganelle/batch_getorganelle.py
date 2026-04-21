#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
batch_getorganelle.py  ——  GetOrganelle 批量脚本（优化版 v7）
Wiki：https://github.com/Kinggerm/GetOrganelle/wiki
FAQ ：https://github.com/Kinggerm/GetOrganelle/wiki/FAQ

v7 优化清单（相较 v6）：
  · [BIOINFO FIX 1] 完成判据收紧：
      - 默认仅将非空 *.path_sequence.fasta 视为可汇总的完整结果
      - graph-only 结果默认归为 INCOMPLETE，避免混入最终 fasta 汇总
      - 新增 --accept_graph_only，允许高级用户显式接受 graph-only 成功状态
  · [BIOINFO FIX 2] INCOMPLETE 独立状态：
      - rc=0 但没有 path_sequence 的样本不再计入 OK
      - 汇总表单独列出 INCOMPLETE，提示用户用 Bandage/日志人工复核
  · [BIOINFO FIX 3] 迭代模式状态修正：
      - --iter_rounds > 1 时跳过判断改为检查最后一轮 iter_N
      - summary 自动只汇总每个样本最高轮次 iter_N，避免同一样本多轮重复入库
  · [BIOINFO FIX 4] 续跑安全 manifest：
      - 每个样本/轮次写入 batch_getorganelle_manifest.json
      - --continue 前检查 reads、seed/label、核心命令是否与上次一致
      - 新增 --ignore_manifest_mismatch 供高级用户强制忽略差异
  · [BIOINFO FIX 5] 失败退出码修正：
      - 默认全流程 assembly 失败后仍可 clean/summary，但最终返回非 0
      - 避免集群任务或上游 workflow 误判为成功
  · [BIOINFO FIX 6] 清理策略收紧：
      - 默认只清理 OK 样本的中间文件
      - FAILED/INCOMPLETE 保留中间文件，便于排错与人工复核
      - 新增 --clean_all_status 允许恢复“无论状态都清理”的行为
  · [BIOINFO FIX 7] 生信运行环境与参数小修：
      - 增加 SPAdes/Bowtie2/BLAST+ 可用性提示
      - -F anonym 缺少 seed/label 时改为错误退出
      - 新增 --subdir_list 正确拼写别名，兼容旧 --sudir_list
      - --force 对已有输出目录追加 GetOrganelle --overwrite
      - GenBank label 自动提取支持 CDS/gene/rRNA/tRNA 输出路径

v6 修复清单（相较 v5）：
  · [BUG FIX 1] --finished_samples 日志语义与统计修正：
      - 日志由"黑名单"改为"已完成样本"，语义准确
      - finished_samples 跳过的样本正确计入汇总表 SKIP 列（v5 中完全丢失）
  · [BUG FIX 2] 子进程卡死修复：
      - subprocess.run(stdout=PIPE) 改为 Popen 逐行流式读取
      - 彻底消除大量输出塞满管道缓冲区导致子进程挂起的问题
  · [BUG FIX 3] 首次运行错误附加 --continue 修复：
      - run_sample() 将 sample_base.mkdir() 移至 resume/done 判断之后
      - is_sample_incomplete() 在脚本自身建目录前执行，
        空目录不再被误判为"已开始未完成"
  · [BUG FIX 4] 参数说明补全：
      - 对照 GetOrganelle 官方文档（Wiki/FAQ）补全所有参数 help 描述
      - 恢复 v4 中完整的 -F 类型列表（embplant_nr / fungus_nr / other_pt）
  · [BUG FIX 5] --no_resume 重新生效：
      - run_sample() 恢复 getattr(args, "resume", True) 检查（v5 中遗漏）
  · [BUG FIX 6] 空任务判断顺序修正：
      - "未找到 reads" 错误判断移至 finished_samples 过滤之前
      - 避免 finished_samples 过滤后 tasks 为空时误报 "未找到 reads"

v5 保留改动（继续沿用）：
  · is_sample_done() 基于真实输出文件判断，不依赖日志内容
  · is_sample_incomplete() 彻底去除"日志有无"判断
  · 完全移除 --overwrite 逻辑（GetOrganelle --continue 对空目录安全）
  · --sample_list / --finished_samples 支持
  · INCOMPLETE graph 不作为 FAILED 处理
"""

import argparse
import json
import logging
import shutil
import subprocess
import sys
import textwrap
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# Logging
# ─────────────────────────────────────────────────────────────────────────────

def setup_logger(log_file: str) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode="a", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )


# ─────────────────────────────────────────────────────────────────────────────
# Subdir list parsing
# ─────────────────────────────────────────────────────────────────────────────

def _read_subdirs_txt(txt: Path) -> list:
    lines = []
    with txt.open(encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                lines.append(s)
    seen, out = set(), []
    for x in lines:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def parse_subdir_list(raw_vals, input_dir: Path):
    if not raw_vals:
        return None
    if len(raw_vals) == 1:
        for candidate in (Path(raw_vals[0]), input_dir / raw_vals[0]):
            if candidate.is_file():
                return _read_subdirs_txt(candidate)
    items = []
    for v in raw_vals:
        for x in v.split(","):
            x = x.strip()
            if x:
                items.append(x)
    seen, out = set(), []
    for x in items:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out or None


# ─────────────────────────────────────────────────────────────────────────────
# sample_list 解析
# ─────────────────────────────────────────────────────────────────────────────

def parse_sample_list(
    txt: Path,
    r1_suffix: str,
    r2_suffix: str,
    input_dir: Path,
) -> "tuple[dict, list]":
    """
    解析 --sample_list 文件。支持三种写法（可混合）：

      写法一：纯样本名
        SampleA
        SampleB
        → 在 input_dir 下查找 SampleA{r1_suffix} / SampleA{r2_suffix}

      写法二：R1 绝对路径
        /data/reads/SampleC_1.clean.fq.gz
        → 自动推导 R2 = /data/reads/SampleC_2.clean.fq.gz
        → 样本名 = 去掉 r1_suffix 的文件名部分

      写法三：混合（上面两种写法可共存）

      # 开头为注释，空行自动忽略

    返回:
      resolved : dict  {sample_name: (fq1_Path, fq2_Path)}  绝对路径条目
      names    : list  [sample_name, ...]                    纯样本名条目
    """
    resolved: dict = {}
    names: list = []
    seen: set = set()

    with txt.open(encoding="utf-8", errors="ignore") as f:
        for lineno, raw_line in enumerate(f, 1):
            s = raw_line.strip()
            if not s or s.startswith("#"):
                continue

            p = Path(s)
            is_path = p.is_absolute() or ("/" in s) or ("\\" in s)

            if is_path:
                fq1 = p.resolve()
                fname = fq1.name

                if not fname.endswith(r1_suffix):
                    logging.warning(
                        f"[SAMPLE_LIST] 行 {lineno}: 文件名不以 R1 后缀 "
                        f"'{r1_suffix}' 结尾，跳过:\n"
                        f"              {fq1}\n"
                        f"              ✎ 请确认 --suffix_fq 参数与文件实际后缀一致"
                    )
                    continue

                sample = fname[: -len(r1_suffix)]
                if not sample:
                    logging.warning(
                        f"[SAMPLE_LIST] 行 {lineno}: 去掉后缀后样本名为空，跳过: {fq1}"
                    )
                    continue

                if sample in seen:
                    logging.warning(
                        f"[SAMPLE_LIST] 行 {lineno}: 重复样本名 '{sample}'，跳过。"
                    )
                    continue

                if not fq1.exists():
                    logging.warning(
                        f"[SAMPLE_LIST] 行 {lineno}: R1 文件不存在，跳过:\n"
                        f"              {fq1}"
                    )
                    continue

                fq2 = fq1.parent / f"{sample}{r2_suffix}"
                if not fq2.exists():
                    logging.warning(
                        f"[SAMPLE_LIST] 行 {lineno}: R2 文件不存在，跳过:\n"
                        f"              {fq2}\n"
                        f"              ✎ 目前只支持 paired-end；单端支持待后续版本"
                    )
                    continue

                seen.add(sample)
                resolved[sample] = (fq1, fq2)
                logging.info(
                    f"[SAMPLE_LIST] 绝对路径 → {sample}:\n"
                    f"              R1: {fq1}\n"
                    f"              R2: {fq2}"
                )

            else:
                sample = s
                if sample in seen:
                    logging.warning(
                        f"[SAMPLE_LIST] 行 {lineno}: 重复样本名 '{sample}'，跳过。"
                    )
                    continue
                seen.add(sample)
                names.append(sample)

    return resolved, names


# ─────────────────────────────────────────────────────────────────────────────
# finished_samples 解析
# ─────────────────────────────────────────────────────────────────────────────

def parse_finished_samples(txt: Path) -> set:
    """
    解析 --finished_samples 文件，返回样本名集合。
    每行一个样本名；# 开头为注释，空行忽略。

    该列表中的样本即使指定 --force 也不会重跑（优先级最高）。
    常用场景：
      · 之前在别处跑完、结果已手动复制过来
      · 质检后锁定不再重跑
      · 配合 --force 保护特定样本
    """
    finished: set = set()
    with txt.open(encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                finished.add(s)
    return finished


# ─────────────────────────────────────────────────────────────────────────────
# Sample detection
# ─────────────────────────────────────────────────────────────────────────────

def detect_samples(flat_dir: Path, r1_suffix: str, r2_suffix: str) -> dict:
    samples = {}
    for fq1 in sorted(flat_dir.glob(f"*{r1_suffix}")):
        if len(fq1.name) <= len(r1_suffix):
            continue
        sample = fq1.name[: -len(r1_suffix)]
        fq2 = flat_dir / f"{sample}{r2_suffix}"
        if fq2.exists():
            samples[sample] = (fq1.resolve(), fq2.resolve())
        else:
            logging.warning(
                f"[SKIP] {flat_dir.name}/{sample}: 未找到 R2 ({fq2.name})"
            )
    return samples


# ─────────────────────────────────────────────────────────────────────────────
# Result dir / fasta discovery
# ─────────────────────────────────────────────────────────────────────────────

_INTER_DIRNAMES = {"seed", "filtered_spades", "extended_spades"}


def _is_intermediate(path: Path) -> bool:
    return any(part in _INTER_DIRNAMES for part in path.parts)


def find_result_dirs(root: Path) -> list:
    found, seen = [], set()

    def _add(d: Path):
        key = str(d.resolve())
        if key not in seen and not _is_intermediate(d):
            seen.add(key)
            found.append(d)

    for p in root.rglob("get_org.log.txt"):
        if p.is_file():
            _add(p.parent)
    if not found:
        for p in root.rglob("*.path_sequence.fasta"):
            if p.is_file() and not _is_intermediate(p.parent):
                _add(p.parent)
    return sorted(found, key=lambda x: str(x))


def list_result_fastas(result_dir: Path, final_only: bool = False) -> list:
    ps_fastas = sorted(result_dir.glob("*.path_sequence.fasta"))
    if ps_fastas:
        return ps_fastas
    # v7/BIOINFO FIX 1：
    # summary/迭代 seed 默认只接受 path_sequence，避免把临时 fasta 或参考库混入。
    if final_only:
        return []
    others = []
    for pat in ("*.fasta", "*.fa", "*.fna"):
        others.extend(sorted(result_dir.glob(pat)))
    seen, uniq = set(), []
    for p in others:
        k = str(p.resolve())
        if k not in seen:
            seen.add(k)
            uniq.append(p)
    return uniq


def _iter_part_index(part: str) -> "int | None":
    if part.startswith("iter_") and part[5:].isdigit():
        return int(part[5:])
    return None


def select_summary_result_dirs(result_dirs: list, output_dir: Path) -> list:
    """
    从所有结果目录中挑选 summary 应使用的目录。

    v7/BIOINFO FIX 3：
      --iter_rounds 产生 sample/iter_0、sample/iter_1... 时，summary 只应取最高
      轮次，否则同一样本会被重复写入 summary_all.fasta。
    """
    non_iter = []
    best_iter = {}

    for d in result_dirs:
        try:
            rel_parts = d.relative_to(output_dir).parts
        except ValueError:
            rel_parts = (d.name,)

        iter_pos, iter_no = None, None
        for i, part in enumerate(rel_parts):
            idx = _iter_part_index(part)
            if idx is not None:
                iter_pos, iter_no = i, idx
                break

        if iter_pos is None:
            non_iter.append(d)
            continue

        sample_key = rel_parts[:iter_pos]
        if not sample_key:
            sample_key = (d.parent.name,)
        old = best_iter.get(sample_key)
        if old is None or iter_no > old[0]:
            best_iter[sample_key] = (iter_no, d)

    selected = []
    iter_keys = set(best_iter)
    for d in non_iter:
        try:
            rel_parts = d.relative_to(output_dir).parts
        except ValueError:
            rel_parts = (d.name,)
        if rel_parts in iter_keys:
            logging.info(
                f"[SUMMARY] 跳过样本根目录（存在 iter_N 结果，取最高轮次）: {d}"
            )
            continue
        selected.append(d)

    for sample_key, (iter_no, d) in sorted(best_iter.items(), key=lambda x: x[0]):
        logging.info(
            f"[SUMMARY] 迭代样本 {'/'.join(sample_key)} 取最高轮次 iter_{iter_no}: {d}"
        )
        selected.append(d)

    return sorted(selected, key=lambda x: str(x))


# ─────────────────────────────────────────────────────────────────────────────
# 判断样本是否已成功完成 / 处于断点状态
#
# GetOrganelle 输出文件层次：
#   get_org.log.txt           —— 日志，运行一开始即创建，不能作为完成判据
#   *.fastg                   —— 组装图（assembly graph），graph 阶段产出
#   *.selected_graph.gfa       —— GFA 格式组装图（新版本 GetOrganelle 常见）
#   *.path_sequence.fasta     —— 完整路径序列，disentangle 成功后才产出（完成标志）
#   seed/ filtered_spades/ extended_spades/  —— 中间目录
#
# v7/BIOINFO FIX 1：
#   默认完成判据收紧为非空 *.path_sequence.fasta。
#   graph-only 结果对植物线粒体/复杂重复样本可能有价值，但不应默认混入
#   “最终序列”汇总，因此仅在 --accept_graph_only 时才视为可跳过完成。
# ─────────────────────────────────────────────────────────────────────────────

_SUCCESS_KEYWORDS = (
    "Finished getting organelle",   # 标准完成行
    "Result status",                # 结果状态行（含 circular / linear）
    "Writing output finished",      # 输出写出完毕
    "Total cost",                   # 总耗时行（每次成功必打印）
    "total cost",
    "succeeded",
    "Organelle assembly finished",
)


def _log_has_success(log_file: Path) -> bool:
    """检查日志是否含成功关键词。"""
    try:
        text = log_file.read_text(encoding="utf-8", errors="ignore")
        return any(kw in text for kw in _SUCCESS_KEYWORDS)
    except Exception:
        return False


def has_path_sequence(run_out: Path) -> bool:
    """是否存在可作为最终序列结果的 path_sequence fasta。"""
    if not run_out.exists():
        return False
    for fa in run_out.glob("*.path_sequence.fasta"):
        if fa.is_file() and fa.stat().st_size > 0:
            return True
    return False


def has_graph_result(run_out: Path) -> bool:
    """是否存在非空组装图。graph 需要人工复核，默认不等同于最终序列。"""
    if not run_out.exists():
        return False
    for pat in ("*.fastg", "*.selected_graph.gfa", "*.gfa"):
        for graph in run_out.glob(pat):
            if graph.is_file() and graph.stat().st_size > 0:
                return True
    return False


def is_graph_only_done(run_out: Path) -> bool:
    """graph-only 且日志显示正常结束；仅供 --accept_graph_only 使用。"""
    log_file = run_out / "get_org.log.txt"
    return has_graph_result(run_out) and log_file.exists() and _log_has_success(log_file)


def is_sample_done(run_out: Path, accept_graph_only: bool = False) -> bool:
    """
    判断样本是否已成功完成。

    v7/BIOINFO FIX 1：
      1. 默认只接受非空 *.path_sequence.fasta；
      2. 若显式 --accept_graph_only，才接受 graph + 成功日志；
      3. 不再仅凭 get_org.log.txt 判断完成，避免把无最终序列样本误判为 OK。
    """
    if has_path_sequence(run_out):
        return True
    if accept_graph_only and is_graph_only_done(run_out):
        return True
    return False


def is_sample_incomplete(run_out: Path, accept_graph_only: bool = False) -> bool:
    """
    判断样本是否处于"已开始但未完成"状态（应使用 --continue 续跑）。

    逻辑：
      - 目录不存在           → False（全新任务，fresh start）
      - 目录存在 + 已完成     → False（无需续跑）
      - 目录存在 + 未完成     → True（传 --continue）

    重要：调用此函数时目录不应由脚本自身预先创建，
          否则空目录会被误判为"已开始未完成"（Bug Fix 3）。
    """
    if not run_out.exists():
        return False
    return not is_sample_done(run_out, accept_graph_only=accept_graph_only)


# ─────────────────────────────────────────────────────────────────────────────
# Clean intermediates
# ─────────────────────────────────────────────────────────────────────────────

_INTER_FILES = [
    "extended_1_paired.fq", "extended_1_unpaired.fq",
    "extended_2_paired.fq", "extended_2_unpaired.fq",
]
_INTER_DIRS = ["filtered_spades", "extended_spades", "seed"]


def clean_intermediates(run_out: Path) -> None:
    for fn in _INTER_FILES:
        fp = run_out / fn
        if fp.is_file():
            try:
                fp.unlink()
                logging.info(f"[CLEAN] 删除: {fp.name}")
            except Exception as e:
                logging.warning(f"[CLEAN] 无法删除 {fp}: {e}")
    for dn in _INTER_DIRS:
        dp = run_out / dn
        if dp.is_dir():
            try:
                shutil.rmtree(dp, ignore_errors=True)
                logging.info(f"[CLEAN] 删除: {dn}/")
            except Exception as e:
                logging.warning(f"[CLEAN] 无法删除 {dp}: {e}")


# ─────────────────────────────────────────────────────────────────────────────
# Runtime environment check
# ─────────────────────────────────────────────────────────────────────────────

def _tool_version_line(tool: str, args=None) -> str:
    args = args or ["--version"]
    try:
        proc = subprocess.run(
            [tool] + args,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=15,
        )
        text = (proc.stdout or "").strip().splitlines()
        return text[0] if text else f"rc={proc.returncode}"
    except Exception as e:
        return f"版本检测失败: {e}"


def check_runtime_environment() -> bool:
    """
    检查 GetOrganelle 及常见依赖是否在 PATH 中。

    v7/BIOINFO FIX 7：
      GetOrganelle 实际依赖 SPAdes/Bowtie2/BLAST+ 等工具。这里不替代官方
      config 检查，但在批量运行前尽早提示环境问题，减少跑到一半才失败。
    """
    required = ["get_organelle_from_reads.py"]
    optional = ["spades.py", "bowtie2", "blastn"]
    version_args = {"blastn": ["-version"]}
    ok = True

    for tool in required:
        path = shutil.which(tool)
        if not path:
            logging.error(f"[ENV] ✗ 未找到必需工具: {tool}")
            ok = False
        else:
            logging.info(f"[ENV] ✓ {tool}: {path}")
            logging.info(
                f"[ENV]   version: {_tool_version_line(tool, version_args.get(tool))}"
            )

    for tool in optional:
        path = shutil.which(tool)
        if not path:
            logging.warning(
                f"[ENV] ? 未在 PATH 中找到 {tool}。"
                "GetOrganelle 运行时可能仍可通过自身配置找到；若失败请检查环境。"
            )
        else:
            logging.info(f"[ENV] ✓ {tool}: {path}")
            logging.info(
                f"[ENV]   version: {_tool_version_line(tool, version_args.get(tool))}"
            )

    return ok


# ─────────────────────────────────────────────────────────────────────────────
# Run manifest for safe resume
# ─────────────────────────────────────────────────────────────────────────────

_MANIFEST_NAME = "batch_getorganelle_manifest.json"


def _path_meta(path: "Path | None"):
    if path is None:
        return None
    try:
        p = Path(path).resolve()
        st = p.stat()
        return {
            "path": str(p),
            "size": st.st_size,
            "mtime": int(st.st_mtime),
        }
    except Exception:
        return {"path": str(path), "missing": True}


def _cmd_for_manifest(cmd: list) -> list:
    # --continue/--overwrite 是运行控制，不代表生物学参数；续跑一致性检查忽略它们。
    return [x for x in cmd if x not in ("--continue", "--overwrite")]


def build_run_manifest(
    group: str,
    sample: str,
    round_index: int,
    fq1: Path,
    fq2: Path,
    run_out: Path,
    cmd: list,
    seed_fasta: "Path | None",
    label_fasta: "Path | None",
    args,
) -> dict:
    return {
        "script": "batch_getorganelle.py",
        "script_version": "v7",
        "group": group,
        "sample": sample,
        "round_index": round_index,
        "run_out": str(run_out.resolve()),
        "command_without_runtime_flags": _cmd_for_manifest(cmd),
        "target_type": str(args.F),
        "inputs": {
            "r1": _path_meta(fq1),
            "r2": _path_meta(fq2),
            "seed": _path_meta(seed_fasta),
            "label": _path_meta(label_fasta),
            "anti_seed": _path_meta(Path(args.anti_seed).resolve()) if getattr(args, "anti_seed", None) else None,
        },
        "core_params": {
            "F": str(args.F),
            "R": str(args.R),
            "t": str(args.t),
            "k": str(args.k),
            "w": getattr(args, "w", None),
            "P": getattr(args, "P", None),
            "max_reads": getattr(args, "max_reads", None),
            "reduce_reads_for_coverage": getattr(args, "reduce_reads_for_coverage", None),
            "max_extending_len": getattr(args, "max_extending_len", None),
            "expected_max_size": getattr(args, "expected_max_size", None),
            "expected_min_size": getattr(args, "expected_min_size", None),
            "min_quality_score": getattr(args, "min_quality_score", None),
            "memory_save": getattr(args, "memory_save", False),
            "all_data": getattr(args, "all_data", False),
            "extra_args": getattr(args, "extra_args", []),
        },
    }


def _manifest_path(run_out: Path) -> Path:
    return run_out / _MANIFEST_NAME


def write_run_manifest(run_out: Path, manifest: dict) -> None:
    try:
        if run_out.exists():
            with _manifest_path(run_out).open("w", encoding="utf-8") as f:
                json.dump(manifest, f, ensure_ascii=False, indent=2, sort_keys=True)
    except Exception as e:
        logging.warning(f"[MANIFEST] 写入失败: {_manifest_path(run_out)}: {e}")


def _manifest_diff(old: dict, new: dict) -> list:
    checks = [
        ("command_without_runtime_flags", old.get("command_without_runtime_flags"), new.get("command_without_runtime_flags")),
        ("inputs", old.get("inputs"), new.get("inputs")),
        ("target_type", old.get("target_type"), new.get("target_type")),
        ("core_params", old.get("core_params"), new.get("core_params")),
    ]
    return [name for name, a, b in checks if a != b]


def manifest_allows_resume(run_out: Path, manifest: dict, args, tag: str, round_label: str) -> bool:
    """
    v7/BIOINFO FIX 4：
      GetOrganelle --continue 不会替用户判断参数是否和上次一致。
      对批量生产而言，reads/seed/label/kmer 等改变后继续旧检查点风险很高。
    """
    mf = _manifest_path(run_out)
    if not mf.exists():
        logging.warning(
            f"[MANIFEST] {tag} {round_label}未找到历史 manifest，"
            "按旧结果目录续跑；建议人工确认参数未变。"
        )
        return True
    try:
        old = json.loads(mf.read_text(encoding="utf-8"))
    except Exception as e:
        if getattr(args, "ignore_manifest_mismatch", False):
            logging.warning(f"[MANIFEST] 读取失败但已忽略: {mf}: {e}")
            return True
        logging.error(f"[MANIFEST] 读取失败，拒绝续跑: {mf}: {e}")
        return False

    diff = _manifest_diff(old, manifest)
    if not diff:
        return True

    if getattr(args, "ignore_manifest_mismatch", False):
        logging.warning(
            f"[MANIFEST] {tag} {round_label}检测到参数/输入变化但已忽略: {', '.join(diff)}"
        )
        return True

    logging.error(
        f"[MANIFEST] {tag} {round_label}检测到参数/输入变化，拒绝 --continue: {', '.join(diff)}\n"
        "           若确认要强行续跑，请加 --ignore_manifest_mismatch；\n"
        "           若要重跑旧输出目录，请使用 --force。"
    )
    return False


# ─────────────────────────────────────────────────────────────────────────────
# PREP_DB helpers
# ─────────────────────────────────────────────────────────────────────────────

def gb_to_fasta_biopython(gb_files: list, out_fasta: Path) -> bool:
    try:
        from Bio import SeqIO
    except ImportError:
        logging.error(
            "[PREP_DB] 未安装 biopython。\n"
            "         请运行：conda install biopython  或  pip install biopython"
        )
        return False
    records = []
    for gb in gb_files:
        try:
            recs = list(SeqIO.parse(str(gb), "genbank"))
            records.extend(recs)
            logging.info(f"[PREP_DB] 读取 {len(recs)} 条记录 from {gb.name}")
        except Exception as e:
            logging.error(f"[PREP_DB] 解析 {gb} 失败: {e}")
            return False
    if not records:
        logging.error("[PREP_DB] GB 文件中未找到任何序列。")
        return False
    try:
        with out_fasta.open("w") as fh:
            SeqIO.write(records, fh, "fasta")
        logging.info(f"[PREP_DB] seed fasta 写入: {out_fasta}（{len(records)} 条）")
        return True
    except Exception as e:
        logging.error(f"[PREP_DB] 写入 {out_fasta} 失败: {e}")
        return False


def cat_fastas(fasta_files: list, out_fasta: Path) -> bool:
    try:
        with out_fasta.open("w") as fout:
            for fp in fasta_files:
                with open(fp, "r", errors="ignore") as fin:
                    content = fin.read()
                    fout.write(content)
                    if content and not content.endswith("\n"):
                        fout.write("\n")
        logging.info(f"[PREP_DB] 合并 {len(fasta_files)} 个 fasta → {out_fasta}")
        return True
    except Exception as e:
        logging.error(f"[PREP_DB] 合并 fasta 失败: {e}")
        return False


def extract_label_from_gb(
    gb_files: list, out_dir: Path, region_type: str = "CDS"
) -> "Path | None":
    out_dir.mkdir(parents=True, exist_ok=True)
    if not shutil.which("get_annotated_regions_from_gb.py"):
        logging.error(
            "[PREP_DB] 未找到 get_annotated_regions_from_gb.py\n"
            "         该工具由 GetOrganelle 提供，请确认已安装并激活对应环境。"
        )
        return None
    cmd = (
        ["get_annotated_regions_from_gb.py"]
        + [str(g) for g in gb_files]
        + ["-o", str(out_dir), "-t", region_type, "--mix", "--overwrite"]
    )
    logging.info(f"[PREP_DB] 提取 label: {' '.join(cmd)}")
    proc = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    if proc.stdout:
        for line in proc.stdout.splitlines():
            logging.info(f"  {line}")
    if proc.returncode != 0:
        logging.error("[PREP_DB] get_annotated_regions_from_gb.py 失败（见上方输出）")
        return None
    # v7/BIOINFO FIX 7：
    # get_annotated_regions_from_gb.py 的输出目录会随 -t 类型变化；
    # 旧版脚本固定查 gene/gene.fasta，导致 -t CDS/rRNA/tRNA 时可能误报失败。
    candidates = [
        out_dir / region_type / f"{region_type}.fasta",
        out_dir / region_type.lower() / f"{region_type.lower()}.fasta",
        out_dir / f"{region_type}.fasta",
        out_dir / f"{region_type.lower()}.fasta",
        out_dir / "gene" / "gene.fasta",
        out_dir / "gene.fasta",
    ]
    label_src = next((p for p in candidates if p.exists()), None)
    if not label_src:
        logging.error(
            "[PREP_DB] 未找到 label fasta，已尝试路径：\n"
            + "\n".join(f"         {p}" for p in candidates)
        )
        return None
    logging.info(f"[PREP_DB] label fasta: {label_src}")
    return label_src


def do_prep_db(args, db_dir: Path):
    db_dir.mkdir(parents=True, exist_ok=True)
    input_files = [Path(f).resolve() for f in args.ref_files]
    if not input_files:
        logging.error("[PREP_DB] 未指定参考文件（--ref_files）。")
        return None, None

    gb_files = [f for f in input_files if f.suffix.lower() in (".gb", ".gbk", ".gbf", ".gbff", ".genbank")]
    fa_files = [f for f in input_files if f.suffix.lower() in (".fasta", ".fa", ".fna", ".fas")]
    unknown  = [f for f in input_files if f not in gb_files and f not in fa_files]
    if unknown:
        logging.warning(f"[PREP_DB] 后缀未识别，当作 fasta: {[x.name for x in unknown]}")
        fa_files.extend(unknown)

    seed_fasta  = db_dir / "seed.fasta"
    label_fasta = None

    if gb_files and fa_files:
        tmp = db_dir / "_tmp_from_gb.fasta"
        if not gb_to_fasta_biopython(gb_files, tmp):
            return None, None
        if not cat_fastas([tmp] + fa_files, seed_fasta):
            return None, None
        tmp.unlink(missing_ok=True)
    elif gb_files:
        if not gb_to_fasta_biopython(gb_files, seed_fasta):
            return None, None
    elif fa_files:
        if not cat_fastas(fa_files, seed_fasta):
            return None, None
    else:
        logging.error("[PREP_DB] 未找到有效的参考文件。")
        return None, None

    logging.info(f"[PREP_DB] seed 数据库: {seed_fasta}")

    if getattr(args, "make_label", False):
        if not gb_files:
            logging.warning("[PREP_DB] --make_label 需要 GB 文件，跳过。")
        else:
            region_type   = getattr(args, "label_region", "CDS")
            label_out_dir = db_dir / "label_regions"
            label_src     = extract_label_from_gb(gb_files, label_out_dir, region_type)
            if label_src:
                label_fasta = db_dir / "label.fasta"
                shutil.copy2(str(label_src), str(label_fasta))
                logging.info(f"[PREP_DB] label 数据库: {label_fasta}")
            else:
                logging.warning("[PREP_DB] label 提取失败。")

    return seed_fasta, label_fasta


# ─────────────────────────────────────────────────────────────────────────────
# --ref_gb auto-preparation
# ─────────────────────────────────────────────────────────────────────────────

def prepare_ref_gb_for_assembly(
    ref_gb_paths: list,
    db_dir: Path,
    region_type: str = "CDS",
    do_seed: bool = True,
    do_label: bool = True,
) -> "tuple[Path | None, Path | None]":
    db_dir.mkdir(parents=True, exist_ok=True)
    if not ref_gb_paths:
        return None, None

    seed_fasta = None
    if do_seed:
        seed_fasta = db_dir / "ref_gb_seed.fasta"
        logging.info(f"[REF_GB] GB → seed fasta: {seed_fasta}")
        if not gb_to_fasta_biopython(ref_gb_paths, seed_fasta):
            return None, None

    label_fasta = None
    if do_label:
        label_out_dir = db_dir / "ref_gb_label_regions"
        label_src     = extract_label_from_gb(ref_gb_paths, label_out_dir, region_type)
        if not label_src:
            logging.warning(
                "[REF_GB] label 提取失败，将只使用 -s seed，不传 --genes。\n"
                "         非 anonym 模式下 GetOrganelle 使用内置 label 库。"
            )
        else:
            label_fasta = db_dir / "ref_gb_label.fasta"
            shutil.copy2(str(label_src), str(label_fasta))
            logging.info(f"[REF_GB] label fasta: {label_fasta}")

    return seed_fasta, label_fasta


# ─────────────────────────────────────────────────────────────────────────────
# Assembly: build command
# ─────────────────────────────────────────────────────────────────────────────

def _build_assembly_cmd(
    fq1: Path,
    fq2: Path,
    run_out: Path,
    args,
    seed_fasta: "Path | None" = None,
    label_fasta: "Path | None" = None,
    resume: bool = False,
    overwrite: bool = False,
) -> list:
    cmd = [
        "get_organelle_from_reads.py",
        "-1", str(fq1),
        "-2", str(fq2),
        "-o", str(run_out),
        "-F", str(args.F),
        "-R", str(args.R),
        "-t", str(args.t),
        "-k", str(args.k),
    ]

    if resume:
        cmd.append("--continue")
    if overwrite:
        # v7/BIOINFO FIX 7：
        # --force 语义是强制重跑；若输出目录已存在，GetOrganelle 需要 --overwrite。
        cmd.append("--overwrite")

    if seed_fasta:
        cmd += ["-s", str(seed_fasta)]
    if label_fasta:
        cmd += ["--genes", str(label_fasta)]

    w_val = getattr(args, "w", None)
    if w_val is not None:
        cmd += ["-w", str(w_val)]

    P_val = getattr(args, "P", None)
    if P_val is not None:
        cmd += ["-P", str(P_val)]

    max_reads = getattr(args, "max_reads", None)
    if max_reads is not None:
        cmd += ["--max-reads", str(max_reads)]

    reduce_cov = getattr(args, "reduce_reads_for_coverage", None)
    if reduce_cov is not None:
        cmd += ["--reduce-reads-for-coverage", str(reduce_cov)]

    mel = getattr(args, "max_extending_len", None)
    if mel is not None:
        cmd += ["--max-extending-len", str(mel)]

    if getattr(args, "memory_save", False):
        cmd.append("--memory-save")

    ems_max = getattr(args, "expected_max_size", None)
    if ems_max is not None:
        cmd += ["--expected-max-size", str(ems_max)]

    ems_min = getattr(args, "expected_min_size", None)
    if ems_min is not None:
        cmd += ["--expected-min-size", str(ems_min)]

    mqs = getattr(args, "min_quality_score", None)
    if mqs is not None:
        cmd += ["--min-quality-score", str(mqs)]

    anti = getattr(args, "anti_seed", None)
    if anti:
        cmd += ["-a", str(anti)]

    prefix = getattr(args, "prefix", None)
    if prefix:
        cmd += ["--prefix", str(prefix)]

    if getattr(args, "all_data", False):
        if not getattr(args, "max_reads", None):
            cmd += ["--max-reads", "inf"]
        if not getattr(args, "reduce_reads_for_coverage", None):
            cmd += ["--reduce-reads-for-coverage", "inf"]

    extra = getattr(args, "extra_args", None)
    if extra:
        cmd.extend(extra)

    return cmd


# ─────────────────────────────────────────────────────────────────────────────
# Assembly: single sample (with iteration)
#
# v6 关键修复：
#   [BUG FIX 2] subprocess.run → Popen 逐行流式读取，防止缓冲区塞满卡死
#   [BUG FIX 3] sample_base.mkdir() 移至 resume/done 判断之后
#   [BUG FIX 5] 恢复 getattr(args, "resume", True) 检查
# ─────────────────────────────────────────────────────────────────────────────

STATUS_OK     = "OK"
STATUS_SKIP   = "SKIP"
STATUS_INCOMPLETE = "INCOMPLETE"
STATUS_FAILED = "FAILED"


def run_sample(
    group: str,
    sample: str,
    fq1: Path,
    fq2: Path,
    args,
    output_dir: Path,
    initial_seed: "Path | None" = None,
    initial_label: "Path | None" = None,
) -> str:
    n_rounds    = getattr(args, "iter_rounds", 1)
    force       = getattr(args, "force", False)
    dry_run     = getattr(args, "dry_run", False)
    accept_graph_only = getattr(args, "accept_graph_only", False)
    tag         = f"{group}/{sample}" if group not in ("", ".", None) else sample
    group_base  = output_dir if group in ("", ".", None) else (output_dir / group)
    sample_base = group_base / sample
    # [BUG FIX 3] 不在此处 mkdir —— 必须在 done/resume 判断之后再建目录

    current_seed  = initial_seed
    current_label = initial_label
    final_status   = STATUS_OK

    for rnd in range(n_rounds):
        run_out     = sample_base if n_rounds == 1 else (sample_base / f"iter_{rnd}")
        round_label = f"轮 {rnd + 1}/{n_rounds} " if n_rounds > 1 else ""

        # ── [BUG FIX 3] done 检查在 mkdir 之前 ──────────────────────────────
        if not force and is_sample_done(run_out, accept_graph_only=accept_graph_only):
            logging.info(
                f"[SKIP] {tag} {round_label}已完成，跳过。"
                "（使用 --force 强制重跑）"
            )
            if rnd < n_rounds - 1:
                fastas = list_result_fastas(run_out, final_only=True)
                if fastas:
                    sample_base.mkdir(parents=True, exist_ok=True)
                    next_seed = sample_base / f"iter_{rnd}_seed.fasta"
                    cat_fastas(fastas, next_seed)
                    current_seed = next_seed
            continue

        # ── [BUG FIX 3 + BUG FIX 5] resume 检查在 mkdir 之前，并恢复 args.resume ──
        # run_out 若不存在 → is_sample_incomplete 返回 False → fresh start（正确）
        # run_out 存在且未完成 → True → 追加 --continue（正确）
        incomplete_flag = is_sample_incomplete(
            run_out, accept_graph_only=accept_graph_only
        )
        resume_flag = (
            getattr(args, "resume", True)
            and incomplete_flag
        )

        if resume_flag:
            logging.info(f"[RESUME] {tag} {round_label}检测到未完成的运行，追加 --continue")
        elif incomplete_flag and not force:
            # v7/BIOINFO FIX 4：
            # --no_resume 的帮助信息声明“由用户手动决定是否清理后重跑”，这里真正执行。
            logging.error(
                f"[FAILED] {tag} {round_label}输出目录已存在但未完成，且 --no_resume 已启用。\n"
                f"         请检查目录: {run_out}\n"
                "         若确认重跑，请使用 --force；若确认续跑，请去掉 --no_resume。"
            )
            return STATUS_FAILED
        else:
            logging.info(f"[RUN] {tag} {round_label}全新任务，直接启动")

        overwrite_flag = force and run_out.exists() and not resume_flag

        # ── 构造命令 ──────────────────────────────────────────────────────────
        cmd = _build_assembly_cmd(
            fq1, fq2, run_out, args,
            current_seed, current_label,
            resume=resume_flag,
            overwrite=overwrite_flag,
        )

        logging.info(f"[CMD] {' '.join(cmd)}")

        if dry_run:
            logging.info("[DRY_RUN] 跳过实际执行（--dry_run 模式）")
            continue

        manifest = build_run_manifest(
            group, sample, rnd, fq1, fq2, run_out, cmd,
            current_seed, current_label, args,
        )

        if resume_flag and not manifest_allows_resume(
            run_out, manifest, args, tag, round_label
        ):
            return STATUS_FAILED

        # ── iter 模式需预建父目录（GetOrganelle 自己创建 run_out）───────────
        if n_rounds > 1:
            sample_base.mkdir(parents=True, exist_ok=True)
        if run_out.exists():
            write_run_manifest(run_out, manifest)

        # ── [BUG FIX 2] Popen 逐行流式读取，避免缓冲区满导致子进程挂起 ──────
        try:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
        except OSError as e:
            logging.error(f"[FAILED] {tag} {round_label}启动失败: {e}")
            return STATUS_FAILED

        for line in proc.stdout:
            logging.info(f"  {line.rstrip()}")
        proc.stdout.close()
        returncode = proc.wait()

        # ── 返回码检查 ────────────────────────────────────────────────────────
        if returncode != 0:
            logging.error(f"[FAILED] {tag} {round_label}(rc={returncode})")
            if run_out.exists():
                write_run_manifest(run_out, manifest)
            return STATUS_FAILED

        if run_out.exists():
            write_run_manifest(run_out, manifest)

        # rc=0 时以文件为准重新检查完成状态
        round_status = STATUS_OK
        if is_sample_done(run_out, accept_graph_only=accept_graph_only):
            logging.info(f"[OK] {tag} {round_label}完成")
        else:
            round_status = STATUS_INCOMPLETE
            final_status = STATUS_INCOMPLETE
            # rc=0 但没有结果文件：incomplete graph（无 path_sequence）
            # v7/BIOINFO FIX 2：
            # GetOrganelle 可能 rc=0 但只给 graph；生信上不能自动等同于最终序列。
            logging.warning(
                f"[INCOMPLETE] {tag} {round_label}"
                f"rc=0 但未找到 *.path_sequence.fasta。\n"
                "  → 组装图可能不完整（incomplete/broken graph），"
                "请用 Bandage 检查 *.fastg。\n"
                "  → 参考：https://github.com/Kinggerm/GetOrganelle/wiki/FAQ"
                "#what-should-i-do-with-incomplete-resultbroken-assembly-graph"
            )

        if (
            getattr(args, "clean_intermediates", True)
            and (round_status == STATUS_OK or getattr(args, "clean_all_status", False))
        ):
            clean_intermediates(run_out)

        if rnd < n_rounds - 1:
            fastas = list_result_fastas(run_out, final_only=True)
            if not fastas:
                logging.warning(f"[ITER] {tag} 轮 {rnd} 无结果 fasta，停止迭代。")
                break
            next_seed = sample_base / f"iter_{rnd}_seed.fasta"
            if not cat_fastas(fastas, next_seed):
                logging.error(f"[ITER] {tag} 合并 fasta 失败，停止迭代。")
                break
            logging.info(f"[ITER] {tag}: 合并 {len(fastas)} 个 fasta → 轮 {rnd+1} seed")
            current_seed = next_seed

    return final_status


# ─────────────────────────────────────────────────────────────────────────────
# Assembly: batch dispatch
# ─────────────────────────────────────────────────────────────────────────────

def do_assembly(
    args,
    input_dir: Path,
    output_dir: Path,
    seed_fasta: "Path | None" = None,
    label_fasta: "Path | None" = None,
) -> int:

    if not getattr(args, "dry_run", False):
        if not check_runtime_environment():
            logging.error(
                "[ASSEMBLY] ✗ 未找到 get_organelle_from_reads.py\n"
                "             请确认：\n"
                "               1. GetOrganelle 已安装（pip install getorganelle）\n"
                "               2. 已激活正确的 conda 环境（conda activate getorganelle）\n"
                "               3. 已初始化数据库（get_organelle_config.py -a embplant_pt）"
            )
            return 1

    r1_suffix, r2_suffix = args.suffix_fq[0], args.suffix_fq[1]
    subdirs = parse_subdir_list(getattr(args, "sudir_list", None), input_dir)

    # ── 解析 sample_list ──────────────────────────────────────────────────────
    sl_path = getattr(args, "sample_list", None)
    sl_resolved: dict = {}
    sl_names:    list = []

    if sl_path:
        sl_txt = Path(sl_path)
        if not sl_txt.is_file():
            logging.error(
                f"[SAMPLE_LIST] ✗ 文件不存在: {sl_txt}\n"
                f"              当前工作目录: {Path.cwd()}"
            )
            return 1
        sl_resolved, sl_names = parse_sample_list(
            sl_txt, r1_suffix, r2_suffix, input_dir
        )
        logging.info(
            f"[SAMPLE_LIST] 解析完成: "
            f"纯样本名 {len(sl_names)} 个，绝对路径 {len(sl_resolved)} 个"
        )

    # ── 解析 finished_samples ─────────────────────────────────────────────────
    fs_path = getattr(args, "finished_samples", None)
    finished_set: set = set()

    if fs_path:
        fs_txt = Path(fs_path)
        if not fs_txt.is_file():
            logging.error(
                f"[FINISHED_SAMPLES] ✗ 文件不存在: {fs_txt}\n"
                f"                   当前工作目录: {Path.cwd()}"
            )
            return 1
        finished_set = parse_finished_samples(fs_txt)
        logging.info(
            f"[FINISHED_SAMPLES] 已加载 {len(finished_set)} 个已完成样本"
            "（即使 --force 也跳过）"
        )

    # ── 是否需要 input_dir ────────────────────────────────────────────────────
    need_input_dir = (sl_path is None) or bool(sl_names) or bool(subdirs)
    if need_input_dir and not input_dir.exists():
        logging.error(
            f"[ASSEMBLY] ✗ input_dir 不存在: {input_dir}\n"
            f"             当前工作目录: {Path.cwd()}"
        )
        return 1

    # ── 构建 tasks ────────────────────────────────────────────────────────────
    tasks: list = []

    if sl_path:
        if sl_names:
            name_set = set(sl_names)
            if subdirs:
                logging.info(f"[ASSEMBLY] 子目录过滤: {subdirs}")
                for sd in subdirs:
                    sd_path = input_dir / sd
                    if not sd_path.is_dir():
                        logging.warning(f"[SKIP] 子目录不存在: {sd_path}")
                        continue
                    for spl, (fq1, fq2) in detect_samples(
                        sd_path, r1_suffix, r2_suffix
                    ).items():
                        if spl in name_set:
                            tasks.append((sd, spl, fq1, fq2))
            else:
                for spl, (fq1, fq2) in detect_samples(
                    input_dir, r1_suffix, r2_suffix
                ).items():
                    if spl in name_set:
                        tasks.append((".", spl, fq1, fq2))

            found_names = {s for _, s, _, _ in tasks}
            for name in sl_names:
                if name not in found_names:
                    logging.warning(
                        f"[SAMPLE_LIST] 样本名 '{name}' 在 input_dir 中未找到对应 reads，"
                        f"跳过。\n"
                        f"              input_dir = {input_dir}\n"
                        f"              R1 后缀   = {r1_suffix}"
                    )

        for spl, (fq1, fq2) in sl_resolved.items():
            tasks.append((".", spl, fq1, fq2))

    else:
        if subdirs:
            logging.info(f"[ASSEMBLY] 子目录过滤: {subdirs}")
            for sd in subdirs:
                sd_path = input_dir / sd
                if not sd_path.is_dir():
                    logging.warning(f"[SKIP] 子目录不存在: {sd_path}")
                    continue
                for spl, (fq1, fq2) in detect_samples(
                    sd_path, r1_suffix, r2_suffix
                ).items():
                    tasks.append((sd, spl, fq1, fq2))
        else:
            for spl, (fq1, fq2) in detect_samples(
                input_dir, r1_suffix, r2_suffix
            ).items():
                tasks.append((".", spl, fq1, fq2))

    # ── [BUG FIX 6] "未找到 reads" 检查在 finished_samples 过滤之前 ──────────
    if not tasks:
        if sl_path:
            logging.error(
                f"[ASSEMBLY] ✗ sample_list 中未找到任何有效的 paired-end 样本\n"
                f"\n"
                f"  当前配置：\n"
                f"    --sample_list = {sl_path}\n"
                f"    --input_dir   = {input_dir}\n"
                f"    R1 后缀       = {r1_suffix}\n"
                f"    R2 后缀       = {r2_suffix}\n"
                f"\n"
                f"  ✎ 解决方法：\n"
                f"    · 样本名写法：检查名称与实际文件名（去掉后缀后）是否一致\n"
                f"    · 绝对路径写法：检查路径是否正确，后缀是否与 --suffix_fq 一致\n"
                f"    · 后缀不符 → 用 --suffix_fq 指定，例如：\n"
                f"        --suffix_fq _R1.fq.gz _R2.fq.gz"
            )
        else:
            try:
                actual_files    = sorted(p.name for p in input_dir.iterdir() if p.is_file())[:30]
                actual_suffixes = sorted({p.suffix for p in input_dir.iterdir() if p.is_file()})
            except Exception:
                actual_files, actual_suffixes = [], []
            logging.error(
                f"[ASSEMBLY] ✗ 未找到 paired-end FASTQ 样本\n"
                f"\n"
                f"  当前配置：\n"
                f"    --input_dir  = {input_dir}\n"
                f"    R1 后缀      = {r1_suffix}\n"
                f"    R2 后缀      = {r2_suffix}\n"
                f"\n"
                f"  目录内实际文件（前 30 个）：\n"
                f"    {actual_files if actual_files else '（目录为空或无法读取）'}\n"
                f"\n"
                f"  目录内文件后缀：{actual_suffixes}\n"
                f"\n"
                f"  ✎ 解决方法：\n"
                f"    · 后缀不符 → 用 --suffix_fq 指定，例如：\n"
                f"        --suffix_fq _R1.fq.gz _R2.fq.gz\n"
                f"        --suffix_fq _1.fastq.gz _2.fastq.gz\n"
                f"    · reads 在子目录 → 用 --subdir_list 指定子目录名\n"
                f"    · 指定部分样本 → 用 --sample_list 指定样本列表文件\n"
                f"    · 路径错误 → 检查 -i 参数"
            )
        return 1

    # ── [BUG FIX 1 + 6] finished_samples 过滤，计入 SKIP ────────────────────
    finished_skipped: list = []
    if finished_set:
        tasks_kept = []
        for g, s, fq1, fq2 in tasks:
            if s in finished_set:
                tag = f"{g}/{s}" if g not in ("", ".", None) else s
                logging.info(f"[DONE] {tag}（已完成样本，跳过）")
                finished_skipped.append(tag)
            else:
                tasks_kept.append((g, s, fq1, fq2))
        if finished_skipped:
            logging.info(
                f"[FINISHED_SAMPLES] 已完成样本共跳过 {len(finished_skipped)} 个"
            )
        tasks = tasks_kept

    n_rounds = getattr(args, "iter_rounds", 1)
    force    = getattr(args, "force", False)
    dry_run  = getattr(args, "dry_run", False)
    accept_graph_only = getattr(args, "accept_graph_only", False)

    # ── 计算各样本状态（跳过 / 待运行）──────────────────────────────────────
    def _group_base(g):
        return output_dir if g in ("", ".", None) else (output_dir / g)

    def _run_out(g, s):
        # v7/BIOINFO FIX 3：
        # 迭代模式是否完成应看最后一轮，而不是 iter_0。
        last_round = max(n_rounds - 1, 0)
        return (
            _group_base(g) / s
            if n_rounds == 1
            else _group_base(g) / s / f"iter_{last_round}"
        )

    tasks_to_run:  list = []
    tasks_skipped: list = []
    for (g, s, fq1, fq2) in tasks:
        ro = _run_out(g, s)
        if not force and is_sample_done(ro, accept_graph_only=accept_graph_only):
            tasks_skipped.append((g, s))
        else:
            tasks_to_run.append((g, s, fq1, fq2))

    total_samples = len(tasks) + len(finished_skipped)
    logging.info(
        f"[ASSEMBLY] 发现 {total_samples} 个样本"
        f"（R1={r1_suffix}  R2={r2_suffix}）\n"
        f"           待运行: {len(tasks_to_run)}  "
        f"已完成跳过: {len(tasks_skipped)}  "
        f"已完成样本跳过: {len(finished_skipped)}  "
        f"并行: {args.max_tasks}  迭代: {n_rounds}"
        + ("  [DRY_RUN]" if dry_run else "")
        + ("  [FORCE]" if force else "")
        + (f"  [SAMPLE_LIST={sl_path}]" if sl_path else "")
        + (f"  [FINISHED_SAMPLES={fs_path}]" if fs_path else "")
    )

    if tasks_skipped:
        for g, s in tasks_skipped:
            tag = f"{g}/{s}" if g not in ("", ".", None) else s
            logging.info(f"[SKIP] {tag}（已完成）")

    if not tasks_to_run:
        logging.info("[ASSEMBLY] 所有样本均已完成，无需重新运行。（使用 --force 强制重跑）")
        results = {
            STATUS_OK: [],
            STATUS_SKIP: [],
            STATUS_INCOMPLETE: [],
            STATUS_FAILED: [],
        }
        for tag in finished_skipped:
            results[STATUS_SKIP].append(tag)
        for g, s in tasks_skipped:
            results[STATUS_SKIP].append(f"{g}/{s}" if g not in ("", ".", None) else s)
        _print_summary_table(results)
        return 0

    # ── 并行执行 ─────────────────────────────────────────────────────────────
    # [BUG FIX 1 + 6] finished_skipped 和 tasks_skipped 都计入 SKIP
    results = {
        STATUS_OK: [],
        STATUS_SKIP: [],
        STATUS_INCOMPLETE: [],
        STATUS_FAILED: [],
    }
    for tag in finished_skipped:
        results[STATUS_SKIP].append(tag)
    for g, s in tasks_skipped:
        results[STATUS_SKIP].append(f"{g}/{s}" if g not in ("", ".", None) else s)

    with ThreadPoolExecutor(max_workers=args.max_tasks) as pool:
        futures = {
            pool.submit(
                run_sample, g, s, fq1, fq2, args,
                output_dir, seed_fasta, label_fasta
            ): (g, s)
            for (g, s, fq1, fq2) in tasks_to_run
        }
        for fut in as_completed(futures):
            g, s = futures[fut]
            try:
                status = fut.result()
            except Exception as exc:
                status = STATUS_FAILED
                logging.error(f"[EXCEPTION] {g}/{s}: {exc}")
            tag = f"{g}/{s}" if g not in ("", ".", None) else s
            logging.info(f"[{status}] {tag}")
            results[status].append(tag)

    _print_summary_table(results)

    return 0 if not results[STATUS_FAILED] and not results[STATUS_INCOMPLETE] else 1


def _print_summary_table(results: dict) -> None:
    ok     = results.get(STATUS_OK,     [])
    skip   = results.get(STATUS_SKIP,   [])
    incomplete = results.get(STATUS_INCOMPLETE, [])
    failed = results.get(STATUS_FAILED, [])
    total  = len(ok) + len(skip) + len(incomplete) + len(failed)
    sep    = "─" * 52

    logging.info("")
    logging.info(sep)
    logging.info(f"  组装汇总（共 {total} 个样本）")
    logging.info(sep)
    logging.info(f"  ✓ 成功   : {len(ok)}")
    logging.info(f"  ⊘ 跳过   : {len(skip)}")
    logging.info(f"  ! 不完整 : {len(incomplete)}")
    logging.info(f"  ✗ 失败   : {len(failed)}")
    if incomplete:
        logging.info("")
        logging.info("  不完整样本列表（需人工复核 graph/log）：")
        for s in incomplete:
            logging.info(f"    · {s}")
    if failed:
        logging.info("")
        logging.info("  失败样本列表：")
        for s in failed:
            logging.info(f"    · {s}")
    logging.info(sep)
    logging.info("")


# ─────────────────────────────────────────────────────────────────────────────
# Clean
# ─────────────────────────────────────────────────────────────────────────────

def do_clean(args, output_dir: Path) -> int:
    if not getattr(args, "clean_intermediates", True):
        logging.info("[CLEAN] --no_clean，跳过。")
        return 0
    result_dirs = find_result_dirs(output_dir)
    logging.info(f"[CLEAN] 找到 {len(result_dirs)} 个结果目录")
    for d in result_dirs:
        clean_intermediates(d)
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

def _read_topology(fa: Path) -> str:
    try:
        with fa.open(encoding="utf-8", errors="ignore") as f:
            for i, line in enumerate(f):
                if line.startswith(">"):
                    low = line.lower()
                    if "circular" in low:
                        return "circular"
                    if "linear" in low or "scaffold" in low:
                        return "linear"
                if i > 500:
                    break
    except Exception:
        pass
    return "linear"


def _write_renamed(src: Path, out_fh, base: str, seq_offset: int, topology: str) -> int:
    seq_idx, last_char = 0, "\n"
    with src.open(encoding="utf-8", errors="ignore") as f:
        for line in f:
            last_char = line[-1] if line else "\n"
            if line.startswith(">"):
                seq_idx += 1
                total_idx = seq_offset + seq_idx
                label = base if total_idx == 1 else f"{base}_{total_idx}"
                out_fh.write(f">{label}  topology={topology}\n")
            else:
                out_fh.write(line)
    if last_char != "\n":
        out_fh.write("\n")
    out_fh.write("\n")
    return seq_idx


def do_summary(output_dir: Path, summary_out: Path) -> int:
    all_result_dirs = find_result_dirs(output_dir)
    result_dirs = select_summary_result_dirs(all_result_dirs, output_dir)
    logging.info(
        f"[SUMMARY] 扫描: {output_dir}  "
        f"结果目录: {len(all_result_dirs)}  汇总目录: {len(result_dirs)}"
    )
    if not result_dirs:
        logging.error("[SUMMARY] 未找到任何结果目录。")
        return 1

    summary_out.parent.mkdir(parents=True, exist_ok=True)
    n_seqs, n_files = 0, 0
    with summary_out.open("w", encoding="utf-8") as out:
        for d in result_dirs:
            try:
                rel_parts = d.relative_to(output_dir).parts
            except ValueError:
                rel_parts = (d.name,)
            sample_parts = [p for p in rel_parts if not p.startswith("iter_")]
            base = "_".join(sample_parts) if sample_parts else d.name

            fastas = list_result_fastas(d, final_only=True)
            if not fastas:
                logging.info(f"[SUMMARY] 跳过（无 path_sequence fasta）: {d}")
                continue
            seq_offset = 0
            for fa in fastas:
                topo    = _read_topology(fa)
                written = _write_renamed(fa, out, base=base, seq_offset=seq_offset, topology=topo)
                seq_offset += written
                n_seqs  += written
                n_files += 1
                logging.info(f"[SUMMARY]  + {base}（{topo}，{written} 条）← {fa.name}")

    if n_seqs == 0:
        # v7/BIOINFO FIX 2：
        # 有日志/graph 但没有 path_sequence 时不能产出有效 summary。
        logging.error("[SUMMARY] 未写入任何 path_sequence 序列，请检查 INCOMPLETE 样本。")
        return 1

    logging.info(f"[SUMMARY] {n_files} 个 fasta / {n_seqs} 条序列 → {summary_out}")
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# Argument parser
# [BUG FIX 4] 补全所有参数说明，恢复 v4 完整 -F 类型列表
# ─────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "-o", "--output_dir", default="getorganelle_results",
        help="批量输出根目录（默认：getorganelle_results）",
    )
    common.add_argument(
        "--log_file", default="getorganelle_batch.log",
        help="日志文件路径（默认：getorganelle_batch.log；追加模式）",
    )

    USAGE = textwrap.dedent("""\
        %(prog)s [子命令] [参数]

        子命令（省略时自动执行 assembly → clean → summary）：
          prep_db   从 GB/FASTA 制备 seed.fasta / label.fasta
          assembly  批量组装（支持迭代、并行、自定义 seed/label）
          clean     清理中间文件
          summary   汇总 *.path_sequence.fasta

        ── 快速开始 ──────────────────────────────────────────────────────────
        %(prog)s -i reads/ -o out/ -F embplant_pt \\
            --suffix_fq _R1.fq.gz _R2.fq.gz

        ── sample_list：只跑指定样本 ─────────────────────────────────────────
        %(prog)s -i reads/ -o out/ -F embplant_pt \\
            --sample_list my_samples.txt

        ── finished_samples：锁定已完成样本（--force 也跳过）────────────────
        %(prog)s -i reads/ -o out/ -F embplant_pt \\
            --finished_samples done.txt

        ── 强制重跑（保护 done.txt 中的样本）───────────────────────────────
        %(prog)s -i reads/ -o out/ -F embplant_pt \\
            --force --finished_samples done.txt

        ── 预览模式：只打印命令，不执行 ─────────────────────────────────────
        %(prog)s -i reads/ -o out/ -F embplant_pt --dry_run

        ── 分步流程（推荐）──────────────────────────────────────────────────
        %(prog)s assembly -i reads/ -o out/ -F embplant_pt --ref_gb ref.gb
        %(prog)s clean    -o out/
        %(prog)s summary  -o out/

        Wiki：https://github.com/Kinggerm/GetOrganelle/wiki
        FAQ ：https://github.com/Kinggerm/GetOrganelle/wiki/FAQ
    """)

    p = argparse.ArgumentParser(
        prog="batch_getorganelle.py",
        usage=USAGE,
        formatter_class=argparse.RawTextHelpFormatter,
        parents=[common],
        add_help=True,
    )

    sub = p.add_subparsers(dest="command", required=False)

    def _add_assembly_args(ap):
        # ── 基础参数 ──────────────────────────────────────────────────────────
        ap.add_argument("-i", "--input_dir", required=True,
                        help="reads 根目录（包含 paired-end FASTQ 文件）")
        ap.add_argument(
            "-F", default="embplant_pt", metavar="TYPE",
            help=textwrap.dedent("""\
                目标基因组类型（默认：embplant_pt）
                  embplant_pt  植物叶绿体（plastome）
                               内置种子库，推荐参数：默认
                  embplant_mt  植物线粒体
                               内置种子库，推荐参数：-R 50
                  embplant_nr  植物核糖体 ITS（18S-ITS1-5.8S-ITS2-28S）
                               内置种子库，推荐参数：-R 10 -k 35,85,105,115
                  animal_mt    动物线粒体
                               内置种子库，推荐参数：默认
                  fungus_mt    真菌线粒体
                               内置种子库，推荐参数：-R 50
                  fungus_nr    真菌核糖体
                               内置种子库，推荐参数：-R 10 -k 35,85,105,115
                  other_pt     其他质体（非植物）
                               内置种子库，推荐参数：默认
                  anonym       完全自定义目标（低拷贝 loci 等）
                               必须同时提供 -s（seed）和 --genes（label）
                               推荐参数：--max_extending_len 100 -P 0"""),
        )

        # ── GetOrganelle 核心参数 ─────────────────────────────────────────────
        go_core = ap.add_argument_group(
            "GetOrganelle 核心参数",
            "对应 get_organelle_from_reads.py 的同名参数，详见官方 Wiki",
        )
        go_core.add_argument(
            "-R", default="15", metavar="INT",
            help=textwrap.dedent("""\
                reads 延伸最大轮次（默认：15）
                  叶绿体通常 10~15 轮即可完成
                  线粒体因基因组更复杂建议设为 50
                  若结果不完整可适当增大"""),
        )
        go_core.add_argument(
            "-t", default="8", metavar="INT",
            help="SPAdes 组装使用的 CPU 线程数（默认：8）",
        )
        go_core.add_argument(
            "-k", default="21,45,65,85,105", metavar="KMERS",
            help=textwrap.dedent("""\
                SPAdes 使用的 kmer 列表，逗号分隔（默认：21,45,65,85,105）
                  核糖体等短区域推荐：35,85,105,115
                  所有 kmer 须为奇数且 ≤ 127"""),
        )
        go_core.add_argument(
            "-w", default=None, metavar="FLOAT_or_INT",
            help=textwrap.dedent("""\
                word size，用于初步 reads 过滤（省略=自动估算）
                  0.0~1.0  → reads 长度的比例（推荐 0.6~0.8）
                  整数 bp  → 固定绝对值
                  过大会导致漏读，过小会引入过多噪声"""),
        )
        go_core.add_argument(
            "-P", default=None, metavar="INT",
            help=textwrap.dedent("""\
                pre-grouping 阶段最多使用的 reads 数（省略=默认 250000）
                  0  → 关闭 pre-grouping（可提高灵敏度，但速度变慢）
                  anonym 模式推荐设为 0"""),
        )
        go_core.add_argument(
            "--max_reads", default=None, metavar="INT_or_inf",
            help=textwrap.dedent("""\
                最多从输入中采样的 reads 数（省略=默认 1500000）
                  inf  → 使用全部 reads（等同于 --all_data 效果之一）
                  对低测序深度样本可设为 inf 以提高组装完整性"""),
        )
        go_core.add_argument(
            "--reduce_reads_for_coverage", default=None, metavar="INT_or_inf",
            help=textwrap.dedent("""\
                达到目标覆盖度后随机削减 reads 的阈值倍数（省略=默认 500）
                  inf  → 关闭削减（等同于 --all_data 效果之一）
                  较小值可防止高覆盖度数据浪费时间"""),
        )
        go_core.add_argument(
            "--max_extending_len", default=None, metavar="INT",
            help=textwrap.dedent("""\
                每轮延伸的最大 reads 长度（省略=不限制）
                  对低拷贝 loci 或 anonym 模式推荐设为 100
                  过大可能引入非目标序列"""),
        )
        go_core.add_argument(
            "--memory_save", action="store_true", default=False,
            help=textwrap.dedent("""\
                启用低内存模式（--memory-save）
                  将 reads 存储到磁盘而非内存，速度略慢
                  适合内存不足（< 8G）的环境"""),
        )
        go_core.add_argument(
            "--expected_max_size", default=None, metavar="INT",
            help=textwrap.dedent("""\
                目标基因组预期最大长度（bp）（--expected-max-size）
                  叶绿体通常 120000~180000
                  超出此值的组装图将被过滤"""),
        )
        go_core.add_argument(
            "--expected_min_size", default=None, metavar="INT",
            help=textwrap.dedent("""\
                目标基因组预期最小长度（bp）（--expected-min-size）
                  低于此值的组装图将被过滤"""),
        )
        go_core.add_argument(
            "--min_quality_score", default=None, metavar="INT",
            help=textwrap.dedent("""\
                reads 碱基质量过滤阈值（Phred 分数，默认：1）
                  建议保持默认；质量过低的碱基将被截断"""),
        )
        go_core.add_argument(
            "--prefix", default=None, metavar="STR",
            help="输出文件名前缀（GetOrganelle --prefix；默认以 -F 类型命名）",
        )
        go_core.add_argument(
            "--anti_seed", dest="anti_seed", default=None, metavar="FASTA",
            help=textwrap.dedent("""\
                反向排除参考序列（-a / --anti-seed）
                  提供与目标基因组高度相似但需排除的序列（如核质 DNA 转移片段）
                  与 seed 格式相同（fasta）"""),
        )

        # ── seed / label 数据库 ───────────────────────────────────────────────
        seed_label = ap.add_argument_group(
            "seed / label 数据库",
            textwrap.dedent("""\
              seed：用于初步 reads 过滤的参考序列（-s）
              label：用于组装图拆解（disentangle）的基因坐标参考（--genes）

              优先级（高 → 低）：
                seed ：-s/--seed > --seed_gb > --ref_gb > 省略（使用内置库）
                label：--genes   >            --ref_gb > 省略（使用内置库）

              注：非 anonym 模式可省略 seed/label，GetOrganelle 使用内置库。"""),
        )
        seed_label.add_argument(
            "-s", "--seed", dest="seed", default=None, metavar="FASTA_or_GB",
            help=textwrap.dedent("""\
                seed 参考序列文件
                  .fasta/.fa/.fna  直接使用
                  .gb/.gbk/.gbf/.gbff/.genbank  自动转换为 fasta
                  可提供近缘物种叶绿体全基因组以提升灵敏度"""),
        )
        seed_label.add_argument(
            "--seed_gb", nargs="+", default=None, metavar="GB",
            help=textwrap.dedent("""\
                从 GenBank 格式文件生成 seed（仅提取序列，不提取 label）
                  若需同时生成 label，请改用 --ref_gb"""),
        )
        seed_label.add_argument(
            "--genes", dest="genes", default=None, metavar="FASTA",
            help=textwrap.dedent("""\
                label 参考序列文件（--genes）
                  包含目标基因组所有基因坐标的 fasta
                  由 get_annotated_regions_from_gb.py 从 GB 文件提取
                  省略时使用内置库（仅 anonym 模式必须提供）"""),
        )
        seed_label.add_argument(
            "--ref_gb", nargs="+", default=None, metavar="GB",
            help=textwrap.dedent("""\
                参考 GenBank 文件，自动生成 seed 和 label（按需）
                  · 若未指定 -s，从 GB 提取全序列作为 seed
                  · 若未指定 --genes，从 GB 提取 CDS/gene 区域作为 label
                  · 同时指定 -s 和 --genes 时，本参数完全忽略"""),
        )
        seed_label.add_argument(
            "--ref_gb_region", default="CDS",
            choices=["CDS", "gene", "rRNA", "tRNA"],
            help="--ref_gb 提取 label 时使用的特征类型（默认：CDS）",
        )

        # ── 运行控制 ──────────────────────────────────────────────────────────
        ctrl = ap.add_argument_group("运行控制")
        ctrl.add_argument(
            "--suffix_fq", nargs=2, metavar=("R1_SUFFIX", "R2_SUFFIX"),
            default=["_1.clean.fq.gz", "_2.clean.fq.gz"],
            help=textwrap.dedent("""\
                reads 文件后缀（默认：_1.clean.fq.gz  _2.clean.fq.gz）
                示例：
                  --suffix_fq _R1.fq.gz _R2.fq.gz
                  --suffix_fq _1.fastq.gz _2.fastq.gz
                  --suffix_fq _R1_001.fastq.gz _R2_001.fastq.gz"""),
        )
        ctrl.add_argument(
            "--max_tasks", type=int, default=4,
            help=textwrap.dedent("""\
                最大并行样本数（默认：4）
                  建议设为 物理核心数 ÷ (-t 线程数)
                  例如 32 核 / 8 线程 = 4 并行"""),
        )
        ctrl.add_argument(
            "--iter_rounds", type=int, default=1,
            help=textwrap.dedent("""\
                迭代组装轮次（默认：1；N>1 时上轮输出作为下轮 seed）
                  适用于组装不完整时通过多轮迭代提升完整性"""),
        )
        ctrl.add_argument(
            "--all_data", action="store_true", default=False,
            help=textwrap.dedent("""\
                使用全量 reads（等价于同时设置：
                  --max_reads inf --reduce_reads_for_coverage inf）
                  适合低测序深度或难组装的样本"""),
        )
        ctrl.add_argument(
            "--subdir_list", "--sudir_list", dest="sudir_list",
            nargs="*", default=None,
            help=textwrap.dedent("""\
                只处理 input_dir 下的指定子目录（空格或逗号分隔，或 txt 文件路径）
                --sudir_list 为历史拼写，仍保留兼容
                示例：
                  --subdir_list batch1 batch2 batch3
                  --subdir_list subdirs.txt"""),
        )
        ctrl.add_argument(
            "--sample_list", "--sl", dest="sample_list",
            default=None, metavar="TXT",
            help=textwrap.dedent("""\
                指定要处理的样本列表文件（每行一个；# 开头为注释，空行忽略）
                支持三种写法（可混合）：
                  写法一：纯样本名（在 --input_dir 下查找 reads）
                    SampleA
                    SampleB
                  写法二：R1 绝对路径（自动推导 R2；可跨目录）
                    /data/project/reads/SampleC_1.clean.fq.gz
                  写法三：混合
                    SampleA
                    /data/other/SampleD_1.clean.fq.gz
                注意：目前只支持 paired-end reads"""),
        )
        ctrl.add_argument(
            "--finished_samples", "-fs", dest="finished_samples",
            default=None, metavar="TXT",
            help=textwrap.dedent("""\
                已完成样本列表文件（每行一个样本名；# 开头为注释，空行忽略）
                该列表中的样本即使指定 --force 也不会重跑（优先级最高）
                常用场景：
                  · 之前在别处跑完、结果已手动复制过来的样本
                  · 质检后锁定不再重跑的样本
                  · 配合 --force 时保护特定样本不被重跑
                示例文件：
                  # 已人工确认的样本
                  SampleA
                  SampleB"""),
        )
        ctrl.add_argument(
            "--no_clean", dest="clean_intermediates",
            action="store_false", default=True,
            help=textwrap.dedent("""\
                保留中间文件（默认会自动清理以节省磁盘空间）
                中间文件包括：
                  filtered_spades/  extended_spades/  seed/
                  extended_*.fq"""),
        )
        ctrl.add_argument(
            "--clean_all_status", action="store_true", default=False,
            help=textwrap.dedent("""\
                无论 OK / INCOMPLETE / FAILED 都清理中间文件
                默认仅清理 OK 样本，保留失败或不完整样本便于排错"""),
        )
        ctrl.add_argument(
            "--no_resume", dest="resume",
            action="store_false", default=True,
            help=textwrap.dedent("""\
                禁用断点续跑（默认对未完成的样本自动追加 --continue）
                指定此选项后：若输出目录已存在但未完成，脚本将直接报错退出
                              由用户手动决定是否清理后重跑"""),
        )
        ctrl.add_argument(
            "--force", action="store_true", default=False,
            help=textwrap.dedent("""\
                强制重跑所有样本（忽略已完成状态；已有输出目录会追加 --overwrite）
                注意：--finished_samples 中的样本不受此参数影响"""),
        )
        ctrl.add_argument(
            "--accept_graph_only", action="store_true", default=False,
            help=textwrap.dedent("""\
                将 graph-only 结果也视为已完成
                默认要求 *.path_sequence.fasta；本选项适合高级用户人工确认 graph 结果后使用"""),
        )
        ctrl.add_argument(
            "--ignore_manifest_mismatch", action="store_true", default=False,
            help=textwrap.dedent("""\
                忽略续跑 manifest 的输入/参数差异
                默认检测到 reads、seed/label、核心参数变化时拒绝 --continue"""),
        )
        ctrl.add_argument(
            "--dry_run", action="store_true", default=False,
            help="预览模式：只打印将要执行的命令，不实际运行",
        )
        ctrl.add_argument(
            "--extra_args", nargs=argparse.REMAINDER, default=[],
            help=textwrap.dedent("""\
                透传给 get_organelle_from_reads.py 的额外参数
                必须放在命令行最后，前面加 -- 分隔
                示例：--extra_args --disentangle-df 10"""),
        )

    # ── prep_db ───────────────────────────────────────────────────────────────
    pd = sub.add_parser(
        "prep_db",
        help="从 GB/FASTA 制备 seed.fasta 和（可选）label.fasta",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=[common],
    )
    pd.add_argument(
        "--ref_files", nargs="+", required=True,
        help=textwrap.dedent("""\
            参考文件列表（支持混合）：
              .gb/.gbk/.gbf/.gbff/.genbank  → 自动转换为 fasta（需 biopython）
              .fasta/.fa/.fna    → 直接合并"""),
    )
    pd.add_argument(
        "--make_label", action="store_true", default=False,
        help="同时从 GB 文件提取基因区域生成 label.fasta（需要 GB 输入）",
    )
    pd.add_argument(
        "--label_region", default="CDS",
        choices=["CDS", "gene", "rRNA", "tRNA"],
        help="提取 label 时使用的特征类型（默认：CDS）",
    )

    # ── assembly ──────────────────────────────────────────────────────────────
    a = sub.add_parser(
        "assembly",
        help="批量运行 GetOrganelle（支持断点续跑、自动跳过已完成样本）",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=[common],
    )
    _add_assembly_args(a)

    # ── clean ─────────────────────────────────────────────────────────────────
    sub.add_parser(
        "clean",
        help="清理所有结果目录中的中间文件",
        parents=[common],
    )

    # ── summary ───────────────────────────────────────────────────────────────
    s = sub.add_parser(
        "summary",
        help="汇总所有 *.path_sequence.fasta 到单个文件",
        parents=[common],
    )
    s.add_argument(
        "--summary_out", default=None,
        help="汇总输出路径（默认：<output_dir>/summary_all.fasta）",
    )

    # ── run ───────────────────────────────────────────────────────────────────
    r = sub.add_parser(
        "run",
        help="显式全流程：assembly → clean → summary",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=[common],
    )
    _add_assembly_args(r)
    r.add_argument(
        "--summary_out", default=None,
        help="汇总输出路径（默认：<output_dir>/summary_all.fasta）",
    )

    # ── 默认全流程 ───────────────────────────────────────────────────────────
    # v7/BIOINFO FIX 7：
    # 旧版把 assembly 的 required -i 同时注册到主解析器，导致
    #   batch_getorganelle.py assembly -i ...
    # 被顶层 required 参数拦截。默认全流程改为 main() 中自动补 "run" 子命令。

    return p


# ─────────────────────────────────────────────────────────────────────────────
# Main helpers
# ─────────────────────────────────────────────────────────────────────────────

def _resolve_summary_out(args, output_dir: Path) -> Path:
    sout = getattr(args, "summary_out", None)
    if not sout:
        return output_dir / "summary_all.fasta"
    p = Path(sout).resolve()
    if p.is_dir() or str(sout).endswith("/"):
        p.mkdir(parents=True, exist_ok=True)
        return p / "summary_all.fasta"
    return p


def _resolve_seed_label(args, output_dir: Path):
    seed_fasta  = Path(args.seed).resolve()  if getattr(args, "seed",  None) else None
    label_fasta = Path(args.genes).resolve() if getattr(args, "genes", None) else None

    if seed_fasta is not None:
        if not seed_fasta.exists():
            logging.error(f"[SEED] ✗ 文件不存在: {seed_fasta}")
            sys.exit(1)
        suffix = seed_fasta.suffix.lower()
        if suffix in (".gb", ".gbk", ".gbf", ".gbff", ".genbank"):
            logging.info(f"[SEED] 检测到 GB 文件（{seed_fasta.name}），自动转换为 fasta …")
            db_dir = output_dir / "db"
            db_dir.mkdir(parents=True, exist_ok=True)
            converted = db_dir / f"{seed_fasta.stem}_seed_converted.fasta"
            if gb_to_fasta_biopython([seed_fasta], converted):
                seed_fasta = converted
            else:
                logging.error(f"[SEED] ✗ GB → fasta 转换失败: {seed_fasta}")
                sys.exit(1)

    seed_gb_files = getattr(args, "seed_gb", None)
    if seed_gb_files:
        if seed_fasta is not None:
            logging.info("[SEED_GB] -s/--seed 已指定，--seed_gb 被忽略。")
        else:
            db_dir    = output_dir / "db"
            db_dir.mkdir(parents=True, exist_ok=True)
            out_fasta = db_dir / "seed_gb.fasta"
            gb_paths  = [Path(g).resolve() for g in seed_gb_files]
            if gb_to_fasta_biopython(gb_paths, out_fasta):
                seed_fasta = out_fasta
            else:
                logging.error("[SEED_GB] GB → fasta 失败，将不传 -s。")

    ref_gb_files = getattr(args, "ref_gb", None)
    if ref_gb_files:
        need_seed  = seed_fasta  is None
        need_label = label_fasta is None
        if not need_seed and not need_label:
            logging.info("[REF_GB] -s 和 --genes 均已指定，--ref_gb 被完全忽略。")
        else:
            db_dir      = output_dir / "db"
            region_type = getattr(args, "ref_gb_region", "CDS")
            gb_paths    = [Path(g).resolve() for g in ref_gb_files]
            auto_seed, auto_label = prepare_ref_gb_for_assembly(
                gb_paths, db_dir, region_type,
                do_seed=need_seed, do_label=need_label,
            )
            if need_seed and auto_seed:
                seed_fasta = auto_seed
            if need_label and auto_label:
                label_fasta = auto_label

    if getattr(args, "F", "") == "anonym":
        missing = []
        if not seed_fasta:
            missing.append("seed（-s / --seed_gb / --ref_gb）")
        if not label_fasta:
            missing.append("label（--genes / --ref_gb）")
        if missing:
            # v7/BIOINFO FIX 7：
            # anonym 模式没有内置目标库可兜底，缺 seed/label 继续运行只会制造无效任务。
            logging.error(
                f"[ASSEMBLY] -F anonym 必须同时提供 {' 和 '.join(missing)}。\n"
                "           低拷贝 loci 还建议：--max_extending_len 100 -P 0"
            )
            sys.exit(1)

    return seed_fasta, label_fasta


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

_SUBCOMMANDS = {"prep_db", "assembly", "clean", "summary", "run"}


def main() -> int:
    parser = build_parser()
    raw = sys.argv[1:]

    if not raw:
        parser.print_help()
        return 0

    _default_mode = (
        raw[0] not in _SUBCOMMANDS
        and raw[0] not in ("-h", "--help")
    )
    parse_raw = ["run"] + raw if _default_mode else raw

    try:
        args = parser.parse_args(parse_raw)
    except SystemExit as e:
        # v7/BIOINFO FIX 7：
        # argparse 在 --help 时会 SystemExit(0)；旧逻辑固定返回 2，导致帮助命令被误判失败。
        return e.code if isinstance(e.code, int) else 2

    setup_logger(args.log_file)
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = args.command if not _default_mode else None

    if cmd == "prep_db":
        db_dir = output_dir / "db"
        seed, label = do_prep_db(args, db_dir)
        return 0 if seed else 1

    if cmd == "assembly":
        input_dir = Path(args.input_dir).resolve()
        seed_fasta, label_fasta = _resolve_seed_label(args, output_dir)
        return do_assembly(args, input_dir, output_dir, seed_fasta, label_fasta)

    if cmd == "clean":
        setattr(args, "clean_intermediates", True)
        return do_clean(args, output_dir)

    if cmd == "summary":
        summary_out = _resolve_summary_out(args, output_dir)
        return do_summary(output_dir, summary_out)

    if cmd in ("run", None):
        if not getattr(args, "input_dir", None):
            parser.print_help()
            return 0
        input_dir   = Path(args.input_dir).resolve()
        summary_out = _resolve_summary_out(args, output_dir)
        seed_fasta, label_fasta = _resolve_seed_label(args, output_dir)

        logging.info("[PIPELINE] Step 1/3: assembly")
        rc = do_assembly(args, input_dir, output_dir, seed_fasta, label_fasta)
        if getattr(args, "dry_run", False):
            # v7/BIOINFO FIX 7：
            # dry-run 只预览 GetOrganelle 命令，不应继续 clean/summary 后因无真实结果失败。
            logging.info("[PIPELINE] --dry_run 模式：跳过 clean + summary。")
            return rc
        if rc != 0:
            logging.warning("[PIPELINE] assembly 有失败样本，继续 clean + summary。")

        logging.info("[PIPELINE] Step 2/3: clean")
        do_clean(args, output_dir)

        logging.info("[PIPELINE] Step 3/3: summary")
        summary_rc = do_summary(output_dir, summary_out)

        # v7/BIOINFO FIX 5：
        # assembly 失败/不完整时，即使 summary 成功也不能让 workflow 误判为成功。
        if rc != 0:
            return rc
        return summary_rc

    logging.error(f"未知子命令: {cmd}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
