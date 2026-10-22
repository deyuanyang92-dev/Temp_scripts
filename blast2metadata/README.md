# blast2metadata

批量序列比对 → 每条 query 取 TopN 命中 → 从 GenBank 抓取命中序列的 metadata（物种、分类、文献…）→ 合并成一张 CSV。
Batch sequence search (BLAST+ local / `blastn -remote` / Biopython qblast / **DIAMOND** / **MMseqs2**)
→ TopN hits per query → GenBank metadata via Entrez → one merged CSV per input file.

典型用途：COI / 16S / ITS 等条形码或线粒体基因的物种鉴定，大批量序列的注释来源追踪。

## 安装 / Install

```bash
# 依赖: Python >= 3.8, Biopython; 可选 pandas, psutil, snakemake
pip install -e /path/to/Temp_scripts/blast2metadata        # 得到 blast2metadata 和 blast2metadata-db 命令
blast2metadata --doctor                                     # 检查 blastn/diamond/mmseqs/snakemake 是否可用

# 不安装也能跑
python /path/to/Temp_scripts/blast2metadata/run_blast2metadata.py -h
```

外部工具按需安装：BLAST+（`blastn`/`blastp`，本地库或 `-remote`）、DIAMOND、MMseqs2、Snakemake。
`conda create -n b2m -c conda-forge -c bioconda blast diamond mmseqs2 biopython pandas psutil snakemake-minimal pytest`

## 快速开始 / Quick start

```bash
# 1) 默认: NCBI 在线 BLAST (blastn -remote) + Entrez metadata
blast2metadata -i query.fasta -o OUT --email you@lab.org

# 2) 本地 nt 库 (BLASTDB 已设置或给路径)
blast2metadata -i query.fasta -o OUT -bm local --blast_db_nt /db/nt --email you@lab.org

# 3) 本地库 + MMseqs2 (核酸或蛋白库; FASTA 会自动 createdb 一次到 OUT/_engine_db/)
blast2metadata -i query.fasta -o OUT -bm local --engine mmseqs --blast_db refs.fasta

# 4) NCBI nr + DIAMOND (蛋白 query -> blastp; 核酸 query -> blastx)。DIAMOND >= 2.1.14 直接读 BLAST 蛋白库，
#    不需要 makedb/.dmnd；库名会像 blastn 一样在 $BLASTDB 里查找
blast2metadata -i proteins.faa -o OUT -bm local --engine diamond --blast_db nr
blast2metadata -i cds.fasta -o OUT -bm local --engine diamond --blast_db /db/nr --query_gencode 5   # 线粒体基因翻译用无脊椎动物密码

# 5) 代理环境下 -remote 常失败 -> 用 Biopython qblast
blast2metadata -i query.fasta -o OUT --online_engine biopython --email you@lab.org

# 6) 只做某一步
blast2metadata -i accessions.tsv -o OUT --mode metadata-only          # accession 列表 -> metadata
blast2metadata -i blast.tsv -o OUT --mode topN-only --topN 5          # 已有 BLAST 表 -> TopN
blast2metadata -i blast.tsv -o OUT --mode add-header --header_style both   # 给 BLAST 表加表头

# 7) 一个目录多个 FASTA (含 .fasta.gz)，默认 Snakemake 调度，可断点续跑
blast2metadata -i fasta_dir -o OUT -bm local --blast_db_nt /db/nt --resume
```

`-h` 显示常用参数，`--help` 显示全部参数。

## 选哪个引擎 / Which engine

| 引擎 | 适用 | 速度 | 注意 |
|---|---|---|---|
| `blast` (默认, BLAST+) | 核酸↔核酸 (blastn/megablast)、蛋白、所有 NCBI 库；`-remote` 在线 | 基线 | 条形码物种鉴定 (高相似度、短序列、对 nt 全库) 仍以 blastn 最准最省事 |
| `diamond` | **只搜蛋白库**：蛋白 query (blastp) 或核酸 query 翻译后 (blastx) | 比 blastp/blastx 快 100–10000× | 不能做 blastn；库可直接是 BLAST 蛋白库 (nr、swissprot，DIAMOND ≥ 2.1.14)、`.dmnd` 或蛋白 FASTA；核酸 query 翻译密码用 `--query_gencode 5` 等 |
| `mmseqs` | 核酸↔核酸、蛋白↔蛋白、核酸→蛋白 (translated) | 大 query 集 (万条以上) 对大库时远快于 BLAST | 每对 query–target 只报一个比对 (不支持核酸多 HSP)，跨环形基因组起点或多外显子的命中 `qcovs` 会偏低；小 query 集每次调用有数秒固定开销 |

DIAMOND/MMseqs2 的输出被整理成与 BLAST+ 相同的 14 列（`qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs`），
`qcovs` 按 BLAST+ 定义由同一 query–subject 所有 HSP 的并集算出，后续 TopN、metadata、表头功能完全一致。
这两个引擎自身多线程，默认 1 个任务占满 `--threads`、`--batch_size 5000`；BLAST+ 默认 4 个并行任务、每批 50 条。

`--blast_db` 可给：BLAST 库名/路径（`blast` 引擎任意库；`diamond` 引擎读蛋白库）；`x.dmnd`；MMseqs2 库前缀（`mmseqs createdb` / `mmseqs databases` 的输出）；或一个 FASTA（自动建库一次到 `OUT/_engine_db/`，`--resume` 和 Snakemake 子任务复用）。

### 下载 nr / nt 直接用 / Getting nr and nt

```bash
# NCBI BLAST 库 (v5 格式, 多卷 nr.00 ... + nr.pal)。BLAST+ 和 DIAMOND >= 2.1.14 都直接用这一份
update_blastdb.pl --decompress --blastdb_version 5 nr      # 蛋白; 同理 nt / core_nt / swissprot
export BLASTDB=/db/blast                                   # 之后 --blast_db nr 即可
blast2metadata -i q.faa -o OUT -bm local --engine diamond --blast_db nr
blast2metadata -i q.fasta -o OUT -bm local --blast_db_nt nt

# MMseqs2 不读 BLAST 库文件，用它自带的下载器 (NR, ClusteredNR, NT, core_nt, UniRef50/90/100, UniProtKB ...)
mmseqs databases NR /db/mmseqs/NR /tmp/mmtmp               # 带 taxonomy
blast2metadata -i q.faa -o OUT -bm local --engine mmseqs --blast_db /db/mmseqs/NR
```

### 用 blast2metadata-db 管理本地库 / Managing local databases

```bash
blast2metadata-db list --remote                  # 列出 $BLASTDB 里的库: 版本/日期/序列数/卷数，并与 NCBI 比较是否有更新
blast2metadata-db list --dir /media/deyuan/18T_A/database/blast_db/nt --dir /media/deyuan/18T_A/database/blast_db/nr --remote
blast2metadata-db update nt -o /db/blast/nt      # 下载 NCBI 清单里的全部 tar.gz (aria2, md5 校验, 断点续传, 失败轮询重试)
                                                  # -> 解压到 .staging_nt/ -> blastdbcmd 校验 -> 替换旧库 -> 删 tar.gz
blast2metadata-db update nr -o /db/blast/nr --downloader ascp        # Aspera (需要 ~/.aspera/connect 和 UDP 33001)
blast2metadata-db update swissprot -o /db/blast --downloader update_blastdb   # 交给 NCBI 的 update_blastdb.pl
blast2metadata-db update nt -o /db/blast/nt --dry_run                # 只打印 URL 列表
blast2metadata-db clean nt -o /db/blast/nt --yes # 删除遗留的 tar.gz / .aria2 / .md5
blast2metadata-db env --dir /db/blast/nt --dir /db/blast/nr          # 打印 export BLASTDB=... 行
```

- 更新判断依据 NCBI 的 `<db>-nucl|prot-metadata.json`（`last-updated`、卷数、大小）对比本地 `blastdbcmd -info` 日期；已最新则不下载（`--force` 强制）。
- 新库先在 `.staging_<db>/` 里解压并通过 `blastdbcmd -info` 校验，才替换旧库（`--keep_old` 把旧库移到 `<db>_old_<日期>/`），中途失败旧库不受影响；重跑同一命令会续传。
- 磁盘空间不足会先报错（需要“下载压缩包 + 解压后”两份空间，`--keep_tarballs` 保留压缩包）。
- 代理：aria2 走 `http_proxy/https_proxy`，`socks://` 的 `all_proxy` 会被忽略；snap 版 aria2c 读不了 `/tmp`，目标目录放 `$HOME` 或 `/media` 下。

DIAMOND 旧文档里的 `diamond prepdb` 已不需要（2.1.14+ 直接 `-d nr`）；`diamond makedb --in nr.gz` 仍可用但没必要。
`blast2metadata --doctor` 会显示你的 DIAMOND 版本是否能直接读 BLAST 库、下载工具 (aria2c/ascp/update_blastdb.pl) 和 `$BLASTDB` 各目录是否存在。

## TopN 排序 / Ranking

- TopN 按 subject 计数（同一 subject 多个 HSP 只留最好的一个），默认按 `bitscore` 排。
- `--sort_by` 预设：`bitscore` `identity` `coverage` `evalue` `pident_x_qcovs`（= pident×qcovs/100，相当于“query 中与 subject 一致的比例”），或自定义列顺序如 `qcovs,pident`。
- `--min_pident` / `--min_qcovs` 先过滤再取 TopN。物种鉴定推荐：`--min_qcovs 80 --sort_by pident_x_qcovs`。
- 不要用 `-max_target_seqs 1` 取“最佳命中”（BLAST 在搜索早期截断，不保证最佳；Shah et al. 2018/2019, Bioinformatics）。

## 输出 / Output

```
OUT/
  <stem>_results.csv         最终表 (UTF-8 BOM, Excel 直接打开)；--header_style code|label|both
  <stem>_completed.txt       完成标记 (--resume 依据)
  blast_results/<stem>_blast.tsv      合并后的比对表 (带表头)；<stem>/batch_*.tsv 各批次
  topN_results/<stem>_topN.tsv
  metadata_results/<stem>_metadata.tsv
  summary_results/<stem>_summary.tsv
  _engine_db/                DIAMOND/MMseqs2 从 FASTA 建的库
  _snakemake_workflow/       Snakefile, config, child_args.json, snakemake_plan.json
  preflight_report.json, batch_blast.log
```

任一批次比对失败 → 该文件失败（非零退出、无完成标记），`--resume` 只重跑失败批次。
Entrez 整批失败（常见原因：某个 accession 无效）会自动拆半重试，其余 accession 仍能拿到 metadata。

## Metadata 来源

- NCBI Entrez（`--email` 必填真实邮箱，`--api_key` 提升到 10 req/s；BLAST 在线接口另有 1 req/3 s 限制，二者独立限速，Snakemake 子任务按并发数均分配额）。
- 本地 TSV：`--local_meta_tsv`，键列自动识别（accession / pseudo_accession / 第一列），`--pseudo_acc_prefixes LOC` 标记本地伪 accession 不去 NCBI 查。
- `--no_ncbi_metadata` 完全离线。

## Snakemake 调度

默认开启（`--no_smk` 关闭）。父进程探测 CPU/内存/负载 (`--smk_busy_guard`) 决定并发数与每任务线程，
每个输入文件一个子任务 `python -m blast2metadata --args_json ... --resume --no_smk`。
相关参数见 `--help` 的 `--smk_*`。

## 开发 / Development

```
src/blast2metadata/
  cli.py        参数解析、main
  pipeline.py   输入发现、预检、单文件流程 (BLAST -> TopN -> metadata)、add-header
  engines.py    BLAST+ 本地/远程、Biopython qblast、DIAMOND、MMseqs2
  metadata.py   本地 metadata TSV、Entrez 抓取、合并
  blasttab.py   BLAST 表解析、合并、TopN 排序
  fields.py     50 个 BLAST 列的字典、表头、-outfmt 处理
  smk.py        Snakemake 计划与 Snakefile
  ratelimit.py  Entrez / BLAST 限速器
  dbtool.py     blast2metadata-db: 本地库清单、与 NCBI 比对、下载更新 (aria2/ascp/update_blastdb.pl)、清理
  accession.py, fastautil.py, writers.py, util.py, doctor.py
tests/          真实 BLAST+/DIAMOND/MMseqs2 回归测试
```

```bash
python -m pytest -q tests                                   # 需要 blastn/makeblastdb 在 PATH
DIAMOND_EXE=/path/diamond MMSEQS_EXE=/path/mmseqs python -m pytest -q tests   # 含引擎测试
B2M_NETWORK_TESTS=1 python -m pytest -q tests -k dbtool      # 真实下载 16S_ribosomal_RNA (~70 MB) 测试更新流程
```

变更记录见 [CHANGELOG.md](CHANGELOG.md)。
