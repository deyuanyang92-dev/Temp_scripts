#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
genome_size.py
批量基因组大小预测：jellyfish count -> jellyfish histo -> GenomeScope2

基于 get_genome_size7.1.py 修改：
  - 改为调用 GenomeScope 2.0（-i/-o/-k/-p 参数），旧版用的是 v1 路径却传了 v2 参数
  - 每步完成写 .done 标记；崩溃后残留的半成品不会被当成"已完成"跳过
  - R1/R2 按样本名配对（旧版 zip 两个未排序的 glob，可能错配）
  - 任一样本失败会在最后汇总报告并以非 0 退出（旧版线程池吞掉异常）
  - 默认跑完 histo 后删除 reads.jf（体积很大），--keep_jf 保留
  - 汇总所有样本的 GenomeScope 结果到 genome_size_summary.tsv

用法：
  # 目录批量（fastp 输出：<dir>/<sample>/<sample>.R1.clean.fastq.gz）
  python3 genome_size.py -i 01.fastp -o 05.genome_size \\
      --genomescope_r /path/to/genomescope2.0/genomescope.R

  # 平铺目录 + 自定义后缀
  python3 genome_size.py -i rawdata -o gs --fq_position flat \\
      --suffix_fq _1.clean.fq.gz _2.clean.fq.gz --genomescope_r ...

  # 单样本
  python3 genome_size.py --fq1 a_1.fq.gz --fq2 a_2.fq.gz --sample a -o gs --genomescope_r ...
"""
import argparse
import os
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def run(cmd, dryrun):
    print(f"[CMD] {cmd}", flush=True)
    if not dryrun:
        subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


def reader_cmd(fq):
    q = shlex.quote(str(fq))
    return f"zcat {q}" if str(fq).endswith(".gz") else f"cat {q}"


def process_pair(sample, fq1, fq2, out_dir, args):
    sd = Path(out_dir) / sample
    sd.mkdir(parents=True, exist_ok=True)
    jf = sd / "reads.jf"
    histo = sd / "reads.histo"
    gs_dir = sd / "genomescope"

    if not (sd / ".histo.done").exists():
        run(f"{args.jellyfish} count -C -m {args.kmer_size} -s {args.hash_size} "
            f"-t {args.threads} -o {shlex.quote(str(jf))} "
            f"<({reader_cmd(fq1)}) <({reader_cmd(fq2)})", args.dryrun)
        run(f"{args.jellyfish} histo -t {args.threads} {shlex.quote(str(jf))} "
            f"> {shlex.quote(str(histo))}", args.dryrun)
        if not args.dryrun:
            (sd / ".histo.done").touch()
            if not args.keep_jf:
                jf.unlink(missing_ok=True)
    else:
        print(f"[SKIP] {sample}: histo 已完成", flush=True)

    if not (sd / ".genomescope.done").exists():
        run(f"{args.rscript} {shlex.quote(args.genomescope_r)} "
            f"-i {shlex.quote(str(histo))} -o {shlex.quote(str(gs_dir))} "
            f"-k {args.kmer_size} -p {args.ploidy} -n {shlex.quote(sample)}",
            args.dryrun)
        if not args.dryrun:
            (sd / ".genomescope.done").touch()
    else:
        print(f"[SKIP] {sample}: GenomeScope 已完成", flush=True)
    return sample


SUMMARY_KEYS = {"Genome Haploid Length": "haploid_length_bp",
                "Heterozygous (ab)": "heterozygosity",
                "Genome Repeat Length": "repeat_length_bp",
                "Genome Unique Length": "unique_length_bp",
                "Model Fit": "model_fit",
                "Read Error Rate": "read_error_rate"}


def find_pairs(reads_dir, fq_position, r1_suffix, r2_suffix):
    """返回 [(sample, fq1, fq2)]，按样本名配对。"""
    root = Path(reads_dir)
    if fq_position == "subdir":
        candidates = [p for d in sorted(root.iterdir()) if d.is_dir()
                      for p in sorted(d.glob(f"*{r1_suffix}"))]
    else:
        candidates = sorted(root.glob(f"*{r1_suffix}"))
    pairs = []
    for r1 in candidates:
        sample = r1.name[: -len(r1_suffix)]
        r2 = r1.with_name(sample + r2_suffix)
        if r2.exists():
            pairs.append((sample, r1, r2))
        else:
            print(f"[WARN] 找不到 R2：{r2}", flush=True)
    return pairs


def parse_summary(gs_dir, sample):
    """读取 GenomeScope2 的 <name>_summary.txt，返回 {col_min: , col_max: }。"""
    f = gs_dir / f"{sample}_summary.txt"
    if not f.exists():
        return None
    row = {}
    for line in f.read_text(errors="ignore").splitlines():
        for k, col in SUMMARY_KEYS.items():
            if line.startswith(k):
                # "1,369,919 bp      1,395,295 bp" -> [min, max]
                vals = [x.replace(",", "") for x in line[len(k):].split() if x != "bp"]
                if len(vals) == 2:
                    row[f"{col}_min"], row[f"{col}_max"] = vals
    return row


def write_summary(out_dir, samples):
    cols = ["sample"] + [f"{c}_{m}" for c in SUMMARY_KEYS.values() for m in ("min", "max")]
    out = Path(out_dir) / "genome_size_summary.tsv"
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for s in samples:
            row = parse_summary(Path(out_dir) / s / "genomescope", s) or {}
            fh.write("\t".join([s] + [row.get(c, "NA") for c in cols[1:]]) + "\n")
    print(f"[OUT] {out}", flush=True)


def main():
    p = argparse.ArgumentParser(description="jellyfish + GenomeScope2 批量基因组大小预测",
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__)
    p.add_argument("-i", "--reads_dir", help="reads 目录")
    p.add_argument("--fq_position", choices=["subdir", "flat"], default="subdir",
                   help="subdir: <dir>/<sample>/reads；flat: <dir>/reads [%(default)s]")
    p.add_argument("--suffix_fq", nargs=2, metavar=("R1", "R2"),
                   default=[".R1.clean.fastq.gz", ".R2.clean.fastq.gz"],
                   help="R1/R2 后缀 [默认匹配 fastp 输出]")
    p.add_argument("--fq1", help="单样本 R1")
    p.add_argument("--fq2", help="单样本 R2")
    p.add_argument("--sample", help="单样本名称")
    p.add_argument("-o", "--out_dir", default="genome_size", help="输出目录")
    p.add_argument("--kmer_size", type=int, default=21, help="k-mer 大小 [%(default)s]")
    p.add_argument("--hash_size", default="20G", help="jellyfish -s [%(default)s]")
    p.add_argument("--threads", type=int, default=16, help="每个样本线程数 [%(default)s]")
    p.add_argument("--max_tasks", type=int, default=2, help="同时运行样本数 [%(default)s]")
    p.add_argument("--ploidy", type=int, default=2, help="GenomeScope2 -p [%(default)s]")
    p.add_argument("--jellyfish", default="jellyfish", help="jellyfish 可执行文件")
    p.add_argument("--rscript", default="Rscript",
                   help="Rscript 可执行文件（需装 genomescope/minpack.lm/argparse 包）")
    p.add_argument("--genomescope_r", required=True, help="GenomeScope2 的 genomescope.R 路径")
    p.add_argument("--keep_jf", action="store_true", help="保留 reads.jf（默认删除）")
    p.add_argument("--dryrun", action="store_true", help="只打印命令")
    args = p.parse_args()

    if args.fq1 and args.fq2:
        sample = args.sample or Path(args.fq1).name.split(".")[0]
        pairs = [(sample, Path(args.fq1), Path(args.fq2))]
    elif args.reads_dir:
        pairs = find_pairs(args.reads_dir, args.fq_position, *args.suffix_fq)
    else:
        p.error("需要 -i/--reads_dir 或 --fq1/--fq2")
    if not pairs:
        sys.exit(f"[ERROR] 未找到 reads，检查 -i / --fq_position / --suffix_fq")
    print(f"[TASKS] {len(pairs)} 个样本: {', '.join(s for s, _, _ in pairs)}", flush=True)

    failed = []
    with ThreadPoolExecutor(max_workers=args.max_tasks) as ex:
        futs = {ex.submit(process_pair, s, r1, r2, args.out_dir, args): s for s, r1, r2 in pairs}
        for fut in as_completed(futs):
            try:
                print(f"[DONE] {fut.result()}", flush=True)
            except Exception as e:
                failed.append(futs[fut])
                print(f"[FAILED] {futs[fut]}: {e}", flush=True)

    if not args.dryrun:
        write_summary(args.out_dir, [s for s, _, _ in pairs])
    if failed:
        sys.exit(f"[ERROR] 失败样本: {', '.join(failed)}")


if __name__ == "__main__":
    main()
