## 2026-10-23 — 在线 BLAST 提速

- Biopython qblast 每批提交一次多条 FASTA，不再逐条提交；按 query ID 核对结果；NCBI 拒绝时对半拆分重试，网络错误不拆分。
- qblast blastn 默认 MEGABLAST=on（与 -remote/本地/NCBI 网页一致），`-megablast off` 可关闭。
- BLAST 提交限速 3.1 s → 10.1 s（NCBI 规定 ≥10 s），`--online_delay_sec` 默认 12 → 0。
- 实测 nt：旧 82.1 s/条；新 50 条一次 72.1 s（1.4 s/条）。数据 `docs/benchmark-online-qblast.jsonl`，脚本 `benchmarks/benchmark_online_qblast.py`。

## 2026-10-23 — 真实 nt 测量后的调整

- 本地 BLAST 默认批大小 500 → 5000：nt 上每批约 13 s 为固定扫库开销，960 条 1 批 66.9 s，20 批 211.9 s，结果相同。
- metadata efetch 每请求 20 → 100 个 accession：200 条实测 4.6–6.4 s 对 17.9–18.5 s。
- TopN 阈值改为作用于每个 subject 的代表（最佳）HSP，短的高一致率片段不再让主比对未达标的 subject 混过 `--min_pident`。
- 新增 `benchmarks/benchmark_nt_passes.py` 与 `docs/benchmark-nt-passes.jsonl`；评估报告补充 nt 实测与 `-mt_mode 1` 不采用的原因。

## 2026-10-22 — BLAST 性能、结果语义与续跑审查

- 修复线程预算向上取整导致超额；本地 BLAST 批次默认 50 → 500，保留显式调参。
- TopN 改为流式、有界候选集合；保留完整 subject ID，稳定处理同分，禁用阈值不要求相应列。
- 严格解析标准 12/14 列、CSV、outfmt 7 和可往返的带引号 TSV；坏行报错，排序仅接受数值列。
- 为搜索、TopN、metadata 和最终结果添加缓存签名；零命中可续跑，输入/参数改变不误用旧结果。
- 原子发布表格/完成标记，外部搜索支持超时，重试清空半成品，qblast 重试撤销已写入的部分 HSP。
- 拒绝重复 query ID 和空序列；修复翻译搜索默认数据库选择和带点本地 metadata ID。
- metadata 不完整时输出缺失清单并保持可续跑状态；Snakemake 将参数文件纳入依赖。
- BLAST 原始表按 QUOTE_NONE 读取（title 以 `"` 开头不再崩溃）；仅读回自写的 TopN TSV 时解析引号。
- 修正 pident×qcovs 的解释；新增回归用例、可运行基准及评估报告，关键改动有中文原因注释。

# Changelog

## 3.0.0 (2026-10-22)

Packaged as `blast2metadata` (`pip install -e .` -> `blast2metadata` command, or
`python -m blast2metadata`, or `python run_blast2metadata.py` from the checkout).
The single 3300-line `blast2_metadatav2.12.py` is split into modules under
`src/blast2metadata/`; the command-line interface is unchanged except for the
additions below. Snakemake child jobs now run `python -m blast2metadata`.

New:
- `--engine diamond|mmseqs` for `-bm local`. DIAMOND (protein DB: blastp /
  blastx) and MMseqs2 (nucleotide or protein DB) write the same 14 columns as
  BLAST+, with `qcovs` computed from the union of HSPs per query-subject as
  BLAST does. `--blast_db` may be a ready `.dmnd` / MMseqs2 DB or a FASTA,
  which is converted once into `OUT/_engine_db/` (lock-protected, reused by
  `--resume` and parallel Snakemake jobs). Batch size defaults to 5000 and one
  job with all threads for these engines.
- `--engine diamond --blast_db nr`: DIAMOND >= 2.1.14 searches a BLAST v5 protein
  database directly (multi-volume `nr.pal` alias included), so the nr from
  `update_blastdb.pl --decompress --blastdb_version 5 nr` needs no makedb/.dmnd;
  bare names are looked up in `$BLASTDB` (DIAMOND itself ignores it), older
  DIAMOND or a nucleotide BLAST DB gives a clear error. MMseqs2 given a BLAST DB
  is pointed at `mmseqs databases NR|NT|core_nt`. `--doctor` reports whether the
  installed DIAMOND can read BLAST DBs.
- `blast2metadata-db` command (list / update / clean / env): lists the databases in
  `$BLASTDB`, compares them with NCBI's `<db>-nucl|prot-metadata.json` manifests,
  downloads a database (aria2c with md5 checksums and retry rounds, ascp, or
  update_blastdb.pl), extracts into a staging directory, verifies with
  `blastdbcmd -info` and swaps it in; cleans leftover tarballs. Replaces the
  ad-hoc aria2_nt.py / ascp_nt.py scripts. `--doctor` now also reports the
  downloaders and whether each `$BLASTDB` directory exists.
- `--query_gencode N` sets the genetic code of translated searches in one place
  (blastx `-query_gencode`, DIAMOND `--query-gencode`, MMseqs2 `--translation-table`).
- gzipped FASTA input (`*.fasta.gz`, output stem drops `.fasta.gz`).
- `--doctor` prints tool/package versions, proxy and BLASTDB; `--version`.
- `--topN`, `--batch_size`, `--max_target_seqs`, `--retries` must be >= 1
  (rejected at parse time).

Robustness:
- A failed Entrez efetch batch (typically one invalid accession -> HTTP 400 for
  the whole batch) is split in halves and retried, so the valid accessions still
  get metadata.
- Ctrl-C cancels the batches not yet started instead of waiting for all of them.
- Biopython qblast now sends the `--email` / `--entrez_tool` (previously only
  Entrez did).
- Duplicate sequence IDs in a query FASTA are reported (their hits pool in TopN).
- Engine executables/DB paths are made absolute before Snakemake children start.

Tests: real-tool tests for DIAMOND and MMseqs2 (skipped when not installed; set
`DIAMOND_EXE` / `MMSEQS_EXE`), gzip input, Entrez batch splitting, entry points.
The outfmt-20 header test is skipped on BLAST+ < 2.13.

## 2.12 and earlier

See the git history of `blast2metadata/blast2_metadatav2.12.py`; the patch
notes that used to head that file (Fix-1 ... Fix-J, Codex 修改, 2026-10-10/11
changes) describe the TopN-by-subject fix, `.part` outputs for safe `--resume`,
the split Entrez/BLAST rate limiters, `--args_json` Snakemake children,
`--mode add-header` and the BLAST field dictionary.
