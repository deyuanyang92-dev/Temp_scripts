#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_pipeline.py
WGS 浅测序统一入口：fastp -> GetOrganelle(线粒体 + 18S-28S) -> SPAdes -> BUSCO -> 基因组大小
                    -> MitoM 线粒体注释 -> BLAST 物种验证（线粒体 + 18S-28S）

本脚本只负责按 config.yaml 依次调用 scripts/ 下的各步骤脚本，各步骤均支持断点续跑
（已完成的样本自动跳过），所以中断后用同一命令重跑即可。

输出目录：
  <out_root>/01.fastp/<sample>/        fastp clean reads + 报告
  <out_root>/02.organelle/             GetOrganelle 线粒体(02_mt) / 18S-28S(03_nr) + 汇总
  <out_root>/03.spades/  04.busco/     SPAdes 组装 / BUSCO 评估
  <out_root>/assembly_busco_summary.tsv
  <out_root>/05.genome_size/           jellyfish + GenomeScope2 + genome_size_summary.tsv
  <out_root>/06.mito_annotation/       MitoM 注释（final_annotation_all.gb, summary.tsv）
  <out_root>/07.blast/                 blast2metadata：TopN 命中 + NCBI 物种信息

用法：
  python3 run_pipeline.py -c config.yaml --dryrun        # 只打印命令
  python3 run_pipeline.py -c config.yaml                 # 全部步骤
  python3 run_pipeline.py -c config.yaml --steps organelle genome_size
"""
import argparse
import shlex
import subprocess
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
SCRIPTS = HERE / "scripts"
STEPS = ["fastp", "organelle", "assembly", "genome_size", "annotate", "blast"]
MT_FASTA, ALL_FASTA = "summary_fasta_mt.fasta", "summary_fasta_all.fasta"
# fastp 输出命名固定为 <sample>/<sample>.R1.clean.fastq.gz
CLEAN_R1, CLEAN_R2 = ".R1.clean.fastq.gz", ".R2.clean.fastq.gz"


def conda_args(cfg, *tools):
    out = []
    for t in tools:
        env = (cfg.get("conda_env") or {}).get(t)
        if env:
            out += [f"--conda_env_{t}", env]
    return out


def reads_args(cfg):
    suffix = str(cfg.get("suffix_fq", "auto")).split()
    return ["-i", cfg["input_dir"], "--fq_position", cfg.get("fq_position", "flat"),
            "--suffix_fq", *suffix]


def cmd_fastp(cfg, out, dry):
    return [sys.executable, SCRIPTS / "fastp_spades_busco.py", "fastp",
            *reads_args(cfg), "-o", out, "--yes",
            "--max_parallel", cfg.get("max_parallel", 4),
            "--threads_fastp", cfg.get("threads", 16),
            *conda_args(cfg, "fastp"), *(["--dryrun"] if dry else []),
            *cfg["extra"].get("fastp", [])]


def cmd_organelle(cfg, out, dry):
    mode = cfg.get("organelle_mode", "all")
    t = cfg.get("threads", 16)
    c = [sys.executable, SCRIPTS / "getorganelle_mt_nr.py", mode,
         "-i", out / "01.fastp", "--fq_position", "subdir",
         "--suffix_fq", f"{CLEAN_R1},{CLEAN_R2}",
         "-o", out / "02.organelle", "--max_parallel", cfg.get("max_parallel", 4),
         "--summary_mt_name", MT_FASTA, "--summary_nr_name", "summary_fasta_nr.fasta",
         "--summary_all_name", ALL_FASTA]
    env = cfg.get("conda_env") or {}
    if env.get("getorganelle"):
        c += ["--conda_env_getorganelle", env["getorganelle"]]
    if mode in ("all", "mgs2get"):
        c += ["--threads_meangs", t, "--species_class_meangs", cfg.get("meangs_species_class", "A-worms")]
        if env.get("meangs"):
            c += ["--conda_env_meangs", env["meangs"]]
    if mode != "nr-only":
        c += ["--t_mt", t]
    if mode in ("all", "mt-nr", "nr-only"):
        c += ["--t_nr", t, "--nr_seed_fasta", cfg["nr_seed_fasta"],
              "--nr_genes_fasta", cfg["nr_genes_fasta"]]
    if dry:
        c.append("--dryrun")
    return c + cfg["extra"].get("organelle", [])


def cmd_assembly(cfg, out, dry):
    # 与 fastp 步骤同一个 out_root：已完成的 fastp 样本会被自动跳过
    return [sys.executable, SCRIPTS / "fastp_spades_busco.py", "all",
            *reads_args(cfg), "-o", out, "--yes",
            "--spades_dir", "03.spades", "--busco_dir", "04.busco",
            "--summary_tsv_name", "assembly_busco_summary.tsv",
            "--busco_lineage", cfg["busco_lineage"],
            "--max_parallel", cfg.get("max_parallel", 4),
            "--threads_fastp", cfg.get("threads", 16),
            "--threads_spades", cfg.get("threads", 16),
            "--threads_busco", cfg.get("threads", 16),
            "--spades_mem_gb", cfg.get("spades_mem_gb", 64),
            *conda_args(cfg, "fastp", "spades", "busco"),
            *(["--dryrun"] if dry else []),
            *cfg["extra"].get("assembly", [])]


def cmd_genome_size(cfg, out, dry):
    return [sys.executable, SCRIPTS / "genome_size.py",
            "-i", out / "01.fastp", "--fq_position", "subdir",
            "--suffix_fq", CLEAN_R1, CLEAN_R2,
            "-o", out / "05.genome_size",
            "--genomescope_r", cfg["genomescope_r"],
            "--rscript", cfg.get("rscript", "Rscript"),
            "--jellyfish", cfg.get("jellyfish", "jellyfish"),
            "--kmer_size", cfg.get("kmer_size", 21),
            "--hash_size", cfg.get("jellyfish_hash", "20G"),
            "--ploidy", cfg.get("ploidy", 2),
            "--threads", cfg.get("threads", 16),
            "--max_tasks", cfg.get("max_parallel", 4),
            *(["--dryrun"] if dry else []),
            *cfg["extra"].get("genome_size", [])]


def cmd_annotate(cfg, out, dry):
    c = [cfg.get("mitom", "mitom"), "run", "-i", out / "02.organelle" / MT_FASTA,
         "-o", out / "06.mito_annotation",
         "--clade", cfg.get("mitom_clade", "auto"),
         "--threads", min(8, int(cfg.get("threads", 16))),
         "--max_tasks", cfg.get("max_parallel", 4)]
    if cfg.get("genetic_code"):
        c += ["--genetic_code", cfg["genetic_code"]]
    if dry:  # mitom 没有 dryrun，只打印命令
        return ["echo", "[DRYRUN]", *c]
    return c + cfg["extra"].get("annotate", [])


def find_blast2metadata(cfg):
    """config 指定 > Temp_scripts 同级目录 > 本目录 external/。"""
    for p in (cfg.get("blast2metadata"),
              HERE.parent / "blast2metadata" / "run_blast2metadata.py",
              HERE / "external" / "blast2metadata" / "run_blast2metadata.py"):
        if p and Path(p).exists():
            return p
    sys.exit("[ERROR] 找不到 run_blast2metadata.py，在 config 中设置 blast2metadata: 路径")


def cmd_blast(cfg, out, dry):
    c = [sys.executable, find_blast2metadata(cfg),
         "-i", out / "02.organelle" / ALL_FASTA, "-o", out / "07.blast",
         "--topN", cfg.get("blast_topN", 5), "--threads", cfg.get("threads", 16),
         "--resume", "--no_smk"]
    if cfg.get("blast_db_nt"):
        c += ["-bm", "local", "--blast_db_nt", cfg["blast_db_nt"]]
    else:  # 在线：biopython qblast 可穿过 HTTP 代理；blastn -remote 在代理下常卡住
        c += ["--online_engine", cfg.get("blast_online_engine", "biopython")]
    if cfg.get("email"):
        c += ["--email", cfg["email"]]
    if dry:
        return ["echo", "[DRYRUN]", *c]
    return c + cfg["extra"].get("blast", [])


BUILDERS = {"fastp": cmd_fastp, "organelle": cmd_organelle,
            "assembly": cmd_assembly, "genome_size": cmd_genome_size,
            "annotate": cmd_annotate, "blast": cmd_blast}
REQUIRED = {"fastp": ["input_dir"], "organelle": [],
            "assembly": ["input_dir", "busco_lineage"], "genome_size": ["genomescope_r"],
            "annotate": [], "blast": ["email"]}


def check_config(cfg, steps):
    missing = [k for s in steps for k in ["out_root"] + REQUIRED[s] if not cfg.get(k)]
    if "organelle" in steps and cfg.get("organelle_mode", "all") in ("all", "mt-nr", "nr-only"):
        missing += [k for k in ("nr_seed_fasta", "nr_genes_fasta") if not cfg.get(k)]
    if missing:
        sys.exit(f"[ERROR] config 缺少: {', '.join(sorted(set(missing)))}")
    for k in ("input_dir", "busco_lineage", "nr_seed_fasta", "nr_genes_fasta", "genomescope_r"):
        if cfg.get(k) and not Path(cfg[k]).exists():
            print(f"[WARN] {k} 路径不存在: {cfg[k]}", flush=True)
    if steps[0] != "fastp" and not (Path(cfg["out_root"]) / "01.fastp").exists() \
            and {"organelle", "genome_size"} & set(steps):
        print("[WARN] 01.fastp 不存在：organelle / genome_size 依赖 fastp clean reads", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-c", "--config", required=True, help="config.yaml")
    p.add_argument("--steps", nargs="+", choices=STEPS, help="覆盖 config 中的 steps")
    p.add_argument("--dryrun", action="store_true", help="各步骤只预览、不运行")
    args = p.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8")) or {}
    cfg["extra"] = cfg.get("extra") or {}
    steps = args.steps or cfg.get("steps") or STEPS
    steps = [s for s in STEPS if s in steps]  # 保证顺序
    check_config(cfg, steps)

    out = Path(cfg["out_root"]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    for step in steps:
        cmd = [str(x) for x in BUILDERS[step](cfg, out, args.dryrun)]
        print(f"\n===== [{step}] =====\n{shlex.join(cmd)}\n", flush=True)
        rc = subprocess.run(cmd).returncode
        # dryrun 时 fastp_spades_busco.py 会因"无输出"返回非 0，忽略
        if rc != 0 and not args.dryrun:
            sys.exit(f"[ERROR] 步骤 {step} 失败 (rc={rc})，修复后用同一命令重跑即可续跑")
    print("\n[DONE] 全部步骤完成。结果在", out, flush=True)


if __name__ == "__main__":
    main()
