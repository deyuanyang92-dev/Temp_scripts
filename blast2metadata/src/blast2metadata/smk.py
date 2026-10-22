"""Snakemake scheduling of per-file jobs."""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Tuple


from .pipeline import process_one_fasta, setup_output_directory
from .util import HostState, _shlex_split, _write_text, detect_host_state, ensure_dir, logger, resolve_total_threads, run_cmd, safe_int, which


def _write_samples_tsv(path: Path, files: List[str]) -> None:
    ensure_dir(path.parent)
    counts: Dict[str, int] = {}
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sample", "fasta"], delimiter="\t")
        w.writeheader()
        for fp in files:
            stem          = Path(fp).stem
            counts[stem]  = counts.get(stem, 0) + 1
            sample        = stem if counts[stem] == 1 else f"{stem}_{counts[stem]}"
            w.writerow({"sample": sample, "fasta": str(Path(fp).resolve())})


def _snakefile_text() -> str:
    """Snakefile text. Child jobs run `python -m blast2metadata --args_json ...`;
    per-job overrides (--input, --output, --out_stem, --resume) are passed explicitly."""
    return r'''
import os, csv, shlex

OUTDIR        = config["outdir"]
PY            = config.get("python_exe", "python")
PKG_ROOT      = config["pkg_root"]
ARGS_JSON     = config["args_json"]
LATENCY_WAIT  = int(config.get("latency_wait", 60))
THREADS_PER_JOB  = int(config["threads_per_job"])
MEM_MB_PER_JOB   = int(config.get("mem_mb_per_job", 0))
JOB_SLOTS_PER_JOB = int(config.get("job_slots_per_job", 1))


def load_samples(tsv):
    rows = []
    with open(tsv, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rows.append(row)
    return rows


SAMPLES    = load_samples(config["samples_tsv"])
BY_SAMPLE  = {r["sample"]: r for r in SAMPLES}
SAMPLE_IDS = [r["sample"] for r in SAMPLES]


rule all:
    input:
        expand(os.path.join(OUTDIR, "{sample}", "{sample}_results.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(OUTDIR, "{sample}", "{sample}_completed.txt"), sample=SAMPLE_IDS)


rule run_one:
    input:
        lambda wc: BY_SAMPLE[wc.sample]["fasta"]
    output:
        csv  = os.path.join(OUTDIR, "{sample}", "{sample}_results.csv"),
        mark = os.path.join(OUTDIR, "{sample}", "{sample}_completed.txt")
    threads:
        THREADS_PER_JOB
    resources:
        mem_mb    = MEM_MB_PER_JOB,
        job_slots = JOB_SLOTS_PER_JOB
    params:
        outdir       = OUTDIR,
        sample_outdir = lambda wc: os.path.join(OUTDIR, wc.sample),
        py           = PY,
        pkg_root     = PKG_ROOT,
        args_json    = ARGS_JSON,
        latency_wait = LATENCY_WAIT
    log:
        os.path.join(OUTDIR, "{sample}", "snakemake.run.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.sample_outdir:q}
        export PYTHONPATH={params.pkg_root:q}${{PYTHONPATH:+:$PYTHONPATH}}
        {params.py:q} -m blast2metadata \
            -i {input:q} \
            -o {params.sample_outdir:q} \
            --args_json {params.args_json:q} \
            --out_stem {wildcards.sample:q} \
            --resume \
            --no_smk \
            > {log:q} 2>&1
        test -s {output.csv:q}
        test -s {output.mark:q}
        """
'''


def _snakemake_error_indicates_lock(stderr: str, stdout: str) -> bool:
    txt     = f"{stderr}\n{stdout}".lower()
    needles = [
        "lockexception", "directory cannot be locked",
        "error: directory cannot be locked",
        "unable to acquire lock", "failed to lock",
    ]
    return any(x in txt for x in needles)


def _run_snakemake_unlock(
    smk: str, snakefile: Path, config_json: Path, wf_dir: Path
) -> Tuple[int, str, str]:
    unlock_cmd = [smk, "-s", str(snakefile), "--configfile", str(config_json), "--unlock"]
    logger.warning("[SMK] Detected lock-related failure, trying automatic --unlock once.")
    return run_cmd(unlock_cmd, cwd=wf_dir, retries=1, retry_sleep=1)


@dataclass
class SnakemakePlan:
    smk_cores:              int
    smk_jobs:               int
    threads_per_job:        int
    total_mem_mb:           int
    mem_mb_per_job:         int
    total_job_slots:        int
    job_slots_per_job:      int
    detected_cpu_count:     int
    detected_mem_total_mb:  int
    detected_mem_available_mb: int
    detected_cpu_percent:   float
    detected_loadavg1:      float
    detected_load_ratio_1m: float
    detected_mem_pressure:  float
    busy_score:             float
    busy_guard_applied:     bool
    busy_cpu_cap:           int
    busy_mem_cap_mb:        int


def _derive_busy_caps(
    host: HostState,
    requested_cores: int,
    requested_mem_mb: int,
    args: argparse.Namespace,
) -> Tuple[int, int, bool]:
    if not args.smk_busy_guard:
        return requested_cores, requested_mem_mb, False

    score = host.busy_score
    if score >= 1.10:
        cpu_frac, avail_mem_frac = 0.15, 0.35
    elif score >= 0.90:
        cpu_frac, avail_mem_frac = 0.25, 0.45
    elif score >= 0.75:
        cpu_frac, avail_mem_frac = 0.35, 0.55
    elif score >= 0.60:
        cpu_frac, avail_mem_frac = 0.50, 0.65
    elif score >= 0.45:
        cpu_frac, avail_mem_frac = 0.65, 0.75
    else:
        cpu_frac      = min(0.90, max(0.30, float(args.smk_idle_cpu_fraction)))
        avail_mem_frac = min(0.90, max(0.35, float(args.smk_idle_avail_mem_fraction)))

    user_cap_avail = float(args.smk_max_available_mem_fraction)
    user_cap_avail = min(0.95, max(0.10, user_cap_avail))
    avail_mem_frac = min(avail_mem_frac, user_cap_avail)

    min_cores  = max(1, safe_int(args.smk_busy_guard_min_cores,  1))
    min_mem_mb = max(1024, safe_int(args.smk_busy_guard_min_mem_mb, 4096))

    cpu_cap    = max(min_cores,  int(host.cpu_count      * cpu_frac))
    mem_cap_mb = max(min_mem_mb, int(host.mem_available_mb * avail_mem_frac))

    final_cores  = max(min_cores,  min(requested_cores,  cpu_cap))
    final_mem_mb = max(min_mem_mb, min(requested_mem_mb, mem_cap_mb))
    busy_applied = final_cores < requested_cores or final_mem_mb < requested_mem_mb
    return final_cores, final_mem_mb, busy_applied


def build_snakemake_plan(args: argparse.Namespace, n_files: int) -> SnakemakePlan:
    host = detect_host_state()
    cpu_count = host.cpu_count

    requested_total_threads = resolve_total_threads(args)
    requested_smk_cores     = safe_int(args.smk_cores, 0)
    smk_cores = (
        max(1, min(requested_smk_cores, cpu_count))
        if requested_smk_cores > 0 else requested_total_threads
    )
    smk_cores = max(1, min(smk_cores, cpu_count))

    total_mem_mb = safe_int(args.smk_total_mem_mb, 0)
    if total_mem_mb <= 0:
        if args.smk_auto_mem:
            base_total = int(
                host.mem_total_mb * min(0.98, max(0.05, float(args.smk_mem_fraction)))
            )
            base_avail = int(
                host.mem_available_mb * min(0.95, max(0.05, float(args.smk_max_available_mem_fraction)))
            )
            total_mem_mb = max(1024, min(base_total, base_avail))
        else:
            total_mem_mb = max(1024, int(host.mem_available_mb * 0.8))

    busy_cpu_cap  = smk_cores
    busy_mem_cap_mb = total_mem_mb
    busy_applied  = False
    if args.smk_busy_guard:
        smk_cores, total_mem_mb, busy_applied = _derive_busy_caps(host, smk_cores, total_mem_mb, args)
        busy_cpu_cap    = smk_cores
        busy_mem_cap_mb = total_mem_mb

    user_jobs           = safe_int(args.smk_jobs, 0)
    user_threads_per_job = safe_int(args.smk_threads_per_job, 0)
    user_mem_per_job    = safe_int(args.smk_mem_mb_per_job, 0)

    auto_min_threads = max(1,    safe_int(args.smk_auto_min_threads_per_job, 4))
    auto_min_mem_mb  = max(1024, safe_int(args.smk_auto_min_mem_mb_per_job, 24576))
    auto_max_jobs    = max(1,    safe_int(args.smk_auto_max_jobs, 8))

    if user_jobs > 0:
        smk_jobs = user_jobs
    else:
        jobs_by_cpu = max(1, smk_cores // auto_min_threads)
        jobs_by_mem = max(1, total_mem_mb // auto_min_mem_mb) if total_mem_mb > 0 else jobs_by_cpu
        smk_jobs    = max(1, min(max(1, n_files), jobs_by_cpu, jobs_by_mem, auto_max_jobs))

    smk_jobs = max(1, min(smk_jobs, max(1, n_files)))

    if user_threads_per_job > 0:
        threads_per_job = user_threads_per_job
    else:
        threads_per_job = max(1, smk_cores // smk_jobs)
        if threads_per_job < auto_min_threads and smk_cores >= auto_min_threads:
            threads_per_job = auto_min_threads

    threads_per_job = max(1, min(threads_per_job, smk_cores))

    max_jobs_by_threads = max(1, smk_cores // threads_per_job)
    if smk_jobs > max_jobs_by_threads:
        smk_jobs = max_jobs_by_threads

    if user_mem_per_job > 0:
        mem_mb_per_job = user_mem_per_job
    elif total_mem_mb > 0:
        mem_mb_per_job = max(1024, int(total_mem_mb / smk_jobs))
    else:
        mem_mb_per_job = 0

    if total_mem_mb > 0 and mem_mb_per_job > 0:
        max_jobs_by_mem = max(1, total_mem_mb // mem_mb_per_job)
        if smk_jobs > max_jobs_by_mem:
            smk_jobs = max_jobs_by_mem
            if user_threads_per_job <= 0:
                threads_per_job = max(1, min(smk_cores, smk_cores // smk_jobs))
            if user_mem_per_job <= 0:
                mem_mb_per_job = max(1024, int(total_mem_mb / smk_jobs))

    total_job_slots   = safe_int(args.smk_total_job_slots, 0)
    if total_job_slots <= 0:
        total_job_slots = smk_jobs
    job_slots_per_job = max(1, safe_int(args.smk_job_slots_per_job, 1))

    return SnakemakePlan(
        smk_cores=smk_cores,
        smk_jobs=smk_jobs,
        threads_per_job=threads_per_job,
        total_mem_mb=total_mem_mb,
        mem_mb_per_job=mem_mb_per_job,
        total_job_slots=total_job_slots,
        job_slots_per_job=job_slots_per_job,
        detected_cpu_count=host.cpu_count,
        detected_mem_total_mb=host.mem_total_mb,
        detected_mem_available_mb=host.mem_available_mb,
        detected_cpu_percent=host.cpu_percent,
        detected_loadavg1=host.loadavg1,
        detected_load_ratio_1m=host.load_ratio_1m,
        detected_mem_pressure=host.mem_pressure,
        busy_score=host.busy_score,
        busy_guard_applied=busy_applied,
        busy_cpu_cap=busy_cpu_cap,
        busy_mem_cap_mb=busy_mem_cap_mb,
    )


def detect_snakemake_cli_capabilities(smk: str) -> Dict[str, object]:
    caps: Dict[str, object] = {"retry_flag": None, "help_ok": False, "version": ""}

    rc, so, se = run_cmd([smk, "--help"], retries=1, retry_sleep=1)
    help_txt = (so or "") + "\n" + (se or "")
    if rc == 0 and help_txt:
        caps["help_ok"] = True
        if "--retries" in help_txt:
            caps["retry_flag"] = "--retries"
        elif "--restart-times" in help_txt:
            caps["retry_flag"] = "--restart-times"
        elif re.search(r"(^|\s)-T([,\s]|$)", help_txt):
            caps["retry_flag"] = "-T"

    vrc, vso, vse = run_cmd([smk, "--version"], retries=1, retry_sleep=1)
    if vrc == 0:
        caps["version"] = (vso or vse or "").strip()
    return caps


def write_smk_plan_json(
    path: Path, plan: SnakemakePlan, caps: Dict[str, object]
) -> None:
    payload = {
        "snakemake_plan": asdict(plan),
        "snakemake_cli":  caps,
        "timestamp":      time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    _write_text(path, json.dumps(payload, indent=2, ensure_ascii=False))


def _args_to_json_safe(args: argparse.Namespace) -> Dict:
    """
    Convert args namespace to a JSON-serialisable dict.
    Non-serialisable values are coerced to str as a fallback.
    """
    d = {}
    for k, v in vars(args).items():
        try:
            json.dumps(v)
            d[k] = v
        except (TypeError, ValueError):
            d[k] = str(v)
    return d


def run_with_snakemake(
    args: argparse.Namespace, files: List[str], out_dir: Path
) -> int:
    smk = which("snakemake")
    if not smk:
        # Codex 修改: --smk_dry_run 在缺少 snakemake 时不再落回真实执行。
        # 理由: 原 fallback 会忽略 dry-run 语义，可能意外启动 BLAST/metadata 请求。
        if args.smk_dry_run:
            logger.error("[SMK] --smk_dry_run requested but snakemake is not installed; not running fallback workflow.")
            return 1
        logger.warning("[SMK] snakemake not found, fallback to python workflow.")
        dirs = setup_output_directory(out_dir)
        # Codex 修改: fallback Python workflow 汇总失败数并返回非零。
        # 理由: 原代码 catch 后仍 return 0，调度系统无法感知失败。
        failed = 0
        for fp in files:
            try:
                if not process_one_fasta(fp, args, dirs):
                    failed += 1
            except Exception as e:
                failed += 1
                logger.error(f"[FILE_FAIL] {fp}: {e}", exc_info=True)
        return 1 if failed else 0

    wf_dir      = ensure_dir(out_dir / "_snakemake_workflow")
    samples_tsv = wf_dir / "samples.tsv"
    config_json = wf_dir / "config.json"
    snakefile   = wf_dir / "Snakefile"
    plan_json   = wf_dir / "snakemake_plan.json"

    # Fix-1: serialise FULL args namespace so child jobs auto-inherit everything
    child_args_json = wf_dir / "child_args.json"
    child_args_dict = _args_to_json_safe(args)
    # Override fields that must differ in child processes
    child_args_dict["use_snakemake"] = False
    # threads_per_job will be overridden below once plan is built; placeholder for now
    _write_text(child_args_json, json.dumps(child_args_dict, indent=2, ensure_ascii=False))

    _write_samples_tsv(samples_tsv, files)
    _write_text(snakefile, _snakefile_text())

    plan = build_snakemake_plan(args, len(files))
    caps = detect_snakemake_cli_capabilities(smk)

    # Update the child args JSON with the actual per-job thread count decided by plan.
    # Fix-A: also inject active-job count so each child process can compute its
    # fair share of the global NCBI rate budget (3 or 10 req/sec ÷ smk_jobs).
    child_args_dict["threads"] = int(plan.threads_per_job)
    child_args_dict["_smk_active_jobs"] = int(plan.smk_jobs)
    _write_text(child_args_json, json.dumps(child_args_dict, indent=2, ensure_ascii=False))

    write_smk_plan_json(plan_json, plan, caps)

    logger.info(
        "[SMK_PLAN] detected_cpu=%s total_mem_mb=%s avail_mem_mb=%s "
        "cpu_percent=%.1f load1=%.2f load_ratio=%.2f mem_pressure=%.2f busy_score=%.2f "
        "busy_guard=%s cores=%s jobs=%s threads_per_job=%s "
        "total_mem_mb=%s mem_mb_per_job=%s total_job_slots=%s job_slots_per_job=%s",
        plan.detected_cpu_count, plan.detected_mem_total_mb, plan.detected_mem_available_mb,
        plan.detected_cpu_percent, plan.detected_loadavg1,
        plan.detected_load_ratio_1m, plan.detected_mem_pressure, plan.busy_score,
        plan.busy_guard_applied, plan.smk_cores, plan.smk_jobs, plan.threads_per_job,
        plan.total_mem_mb, plan.mem_mb_per_job, plan.total_job_slots, plan.job_slots_per_job,
    )
    logger.info("[SMK_CLI] version=%s retry_flag=%s", caps.get("version", ""), caps.get("retry_flag"))

    # Fix-5: config now carries args_json path instead of a long base_args list
    cfg: Dict[str, object] = {
        "outdir":           str(out_dir.resolve()),
        "python_exe":       sys.executable,
        # 修改 (2026-10-22): 子任务用 `python -m blast2metadata`；pkg_root 是包所在目录（src/ 或 site-packages）。
        "pkg_root":         str(Path(__file__).resolve().parent.parent),
        "samples_tsv":      str(samples_tsv.resolve()),
        "args_json":        str(child_args_json.resolve()),
        "threads_per_job":  int(plan.threads_per_job),
        "mem_mb_per_job":   int(plan.mem_mb_per_job),
        "job_slots_per_job": int(plan.job_slots_per_job),
        "latency_wait":     int(args.smk_latency_wait),
    }
    _write_text(config_json, json.dumps(cfg, indent=2, ensure_ascii=False))

    # Fix-B: In Snakemake v8+ `--cores` / `-c` and `--jobs` / `-j` are strict
    # aliases.  Passing both causes: "error: not allowed with argument --jobs/-j".
    # We only pass `--cores`; concurrency is controlled via --resources job_slots.
    cmd = [
        smk,
        "-s", str(snakefile),
        "--configfile", str(config_json),
        "--cores", str(plan.smk_cores),
        "--latency-wait", str(int(args.smk_latency_wait)),
    ]

    retry_flag = str(caps.get("retry_flag") or "").strip()
    if retry_flag:
        cmd += [retry_flag, str(int(args.smk_restart_times))]
    else:
        logger.warning(
            "[SMK] Could not detect a supported retry flag from snakemake --help; "
            "skipping CLI retry option."
        )

    if args.smk_rerun_incomplete:
        cmd.append("--rerun-incomplete")
    if args.smk_keep_going:
        cmd.append("--keep-going")
    if args.smk_printshellcmds:
        cmd.append("--printshellcmds")
    if args.smk_dry_run:
        cmd.append("--dry-run")

    resources_chunks: List[str] = []
    if plan.total_mem_mb > 0 and plan.mem_mb_per_job > 0:
        resources_chunks.append(f"mem_mb={plan.total_mem_mb}")
    if plan.total_job_slots > 0 and plan.job_slots_per_job > 0:
        resources_chunks.append(f"job_slots={plan.total_job_slots}")
    if resources_chunks:
        cmd += ["--resources"] + resources_chunks

    if args.smk_profile:
        cmd += ["--profile", args.smk_profile]
    if args.smk_extra:
        cmd += _shlex_split(args.smk_extra)

    rc, so, se = run_cmd(cmd, cwd=wf_dir, retries=1, retry_sleep=1)
    if rc != 0 and args.smk_unlock_on_lock_exception and _snakemake_error_indicates_lock(se, so):
        urc, uso, use = _run_snakemake_unlock(smk, snakefile, config_json, wf_dir)
        if urc == 0:
            logger.info("[SMK] unlock success, rerunning workflow.")
            rc, so, se = run_cmd(cmd, cwd=wf_dir, retries=1, retry_sleep=1)
        else:
            logger.error(f"[SMK_UNLOCK_FAIL] rc={urc}\n{use[-2000:]}")
            return urc

    if rc != 0:
        logger.error(
            f"[SMK_FAIL] rc={rc}\nSTDOUT_tail:\n{so[-3000:]}\nSTDERR_tail:\n{se[-3000:]}"
        )
        return rc

    logger.info("[SMK] Completed.")
    return 0
