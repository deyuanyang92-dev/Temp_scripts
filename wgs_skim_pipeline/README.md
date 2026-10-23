# WGS 浅测序流程（wgs_skim_pipeline）

每个工具都能**单独使用**；`run_pipeline.py` 按 config.yaml 把它们串成一键流程。

```
原始 reads
  └─ fastp 质控 ............................ 01.fastp/
       ├─ MEANGS 种子 → GetOrganelle 线粒体 + 18S-28S ... 02.organelle/
       │     ├─ MitoM 线粒体注释 ................ 06.mito_annotation/
       │     └─ BLAST 物种验证（mt + nr） ........ 07.blast/
       ├─ SPAdes 组装 → BUSCO 完整度 ........... 03.spades/ 04.busco/
       └─ jellyfish + GenomeScope2 基因组大小 ... 05.genome_size/
```

## 目录

```
wgs_skim_pipeline/
  run_pipeline.py            一键入口（读 config.yaml）
  config.example.yaml        配置模板
  scripts/
    fastp_spades_busco.py    fastp / SPAdes / BUSCO，可单独或组合使用（v6.25 + BUSCO 修复）
    getorganelle_mt_nr.py    MEANGS → GetOrganelle 线粒体 → 18S-28S 一键（v4）
    genome_size.py           jellyfish + GenomeScope2（由 get_genome_size7.1.py 修复改写）
```

同仓库其他目录：`../getorganelle/batch_getorganelle.py`（通用 GetOrganelle：任意 -F、迭代组装、prep_db、summary），
`../blast2metadata/`（BLAST → TopN → NCBI 物种信息，BLAST 步骤自动调用）。

线粒体注释用 MitoM（<https://github.com/deyuanyang92-dev/Mitom>），需单独安装：
`git clone ... && bash install.sh`，得到 `mitom` 命令。

## 一键运行

```bash
cp config.example.yaml config.yaml          # 修改 input_dir / out_root / 数据库路径 / email
python3 run_pipeline.py -c config.yaml --dryrun                 # 预览
python3 run_pipeline.py -c config.yaml                          # 全部步骤
python3 run_pipeline.py -c config.yaml --steps organelle annotate blast   # 部分步骤
```

步骤：`fastp organelle assembly genome_size annotate blast`。中断后用同一命令重跑，已完成样本自动跳过。
建议在 `tmux` / `nohup` 中运行。

## 各工具单独使用

以下 `S=scripts`。所有脚本 `-h` 查看完整参数；`--dryrun` 预览命令。
conda 环境用 `--conda_env_<工具> 环境名` 指定，不写则用 PATH 中的程序。

### 1. fastp（WGS / RNA）

```bash
python3 $S/fastp_spades_busco.py fastp -i rawdata -o out --data_type wgs --conda_env_fastp fastp --yes
python3 $S/fastp_spades_busco.py fastp -i rawdata -o out --data_type rna --conda_env_fastp fastp --yes   # 额外 --trim_poly_x
```

输出 `out/01.fastp/<样本>/<样本>.R1.clean.fastq.gz`、`fastp.html`。

### 2. SPAdes（多种数据类型）

```bash
python3 $S/fastp_spades_busco.py spades -i cleandata -o out --conda_env_spades spades --yes                       # 基因组 (WGS)
python3 $S/fastp_spades_busco.py spades -i cleandata -o out --data_type rna --conda_env_spades spades --yes      # 转录组 --rna
python3 $S/fastp_spades_busco.py spades -i cleandata -o out --assembly_mode meta --conda_env_spades spades --yes  # 宏基因组
```

`--assembly_mode`：isolate sc meta plasmid metaplasmid metaviral rna rnaviral bio corona sewage。
中断后重跑自动 `spades.py --continue`；默认完成后清理中间文件（`--no-spades_clean` 保留）。

### 3. BUSCO（组装完整度）

```bash
# 对已有组装（目录下每个样本一个 fasta）
python3 $S/fastp_spades_busco.py busco -o out --busco_fasta_dir assemblies/ \
    --busco_lineage /media/deyuan/18T_A/database/busco/metazoa_odb12 --conda_env_busco busco --yes
# SPAdes → BUSCO
python3 $S/fastp_spades_busco.py spades2busco -i cleandata -o out --busco_lineage ... --yes
# fastp → SPAdes → BUSCO
python3 $S/fastp_spades_busco.py all -i rawdata -o out --busco_lineage ... --yes
```

RNA 数据自动用 `-m transcriptome`。汇总 `out/summary_report.tsv`。其他：`status`（查看进度）、`clean`、`check`（检查软件）。

### 4. GetOrganelle：MEANGS 种子 + 线粒体 + 18S-28S 一键

```bash
python3 $S/getorganelle_mt_nr.py all -i 01.fastp --fq_position subdir \
    --suffix_fq ".R1.clean.fastq.gz,.R2.clean.fastq.gz" -o organelle \
    --nr_seed_fasta /media/deyuan/18T_A/database/18-28s/seed.fasta \
    --nr_genes_fasta /media/deyuan/18T_A/database/18-28s/lable.fasta
```

| 子命令 | 内容 |
|---|---|
| `all` | MEANGS → 线粒体（MEANGS 结果作 `-s`）→ 18S-28S |
| `mgs2get` | MEANGS → 线粒体 |
| `mt-nr` | 线粒体（GetOrganelle 内置种子）+ 18S-28S |
| `getorganelle` | 只线粒体 |
| `nr-only` | 只 18S-28S（`-F anonym -s seed --genes label -P 0`） |

输出 `summary_fasta_mt/nr/all`、`summary_report.tsv`（长度、环状、N 数）、`failed_all.txt`。

### 5. GetOrganelle 通用 / 迭代组装

```bash
python3 ../getorganelle/batch_getorganelle.py prep_db --ref_files ref.gb --make_label -o db/          # GenBank → seed/label
python3 ../getorganelle/batch_getorganelle.py assembly -i reads/ -o out/ -F animal_mt -s seed.fasta --iter_rounds 3
python3 ../getorganelle/batch_getorganelle.py summary -o out/
```

`--iter_rounds N`：上一轮结果作下一轮种子。适合种子与目标物种差异大的情况。

### 6. 基因组大小

```bash
python3 $S/genome_size.py -i out/01.fastp -o genome_size \
    --genomescope_r /home/deyuan/Softs/genomescope2.0-master/genomescope.R \
    --rscript /home/deyuan/miniforge3/envs/genome_size/bin/Rscript
```

输出 `genome_size_summary.tsv`（基因组单倍体长度、杂合度、重复长度、拟合度的 min/max）和每样本 GenomeScope 图。

### 7. 线粒体注释（MitoM）

```bash
mitom run -i organelle/summary_fasta_mt.fasta -o mito_anno --clade Annelida-segmented-worms --genetic_code 5
```

### 8. BLAST 物种 / 污染验证

```bash
# 在线（默认，NCBI qblast，无需本地库）
python3 ../blast2metadata/run_blast2metadata.py -i organelle/summary_fasta_all.fasta -o blast \
    --online_engine biopython --topN 5 --email you@lab.org --no_smk
# 本地 nt（样本多时更快；NCBI 在线限速约 10 秒/次请求）
python3 ../blast2metadata/run_blast2metadata.py -i organelle/summary_fasta_all.fasta -o blast \
    -bm local --blast_db_nt /media/deyuan/18T_A/database/blast_db/nt/nt --topN 5 --email you@lab.org --no_smk
```

一键流程中由 config 的 `blast_db_nt` 控制：`null` = 在线，填库前缀 = 本地。
在线默认用 Biopython qblast（`blast_online_engine`）：本机有 HTTP 代理，`blastn -remote` 测试时卡住 25 分钟无结果。
本地全 nt 对 16 kb 线粒体序列较慢（测试中 20 分钟未完成），样本少时在线更省事。

查看每条序列 TopN 命中的物种：与预期物种不符 → 鉴定有误或污染。

## 依赖（本机 conda 环境）

| 工具 | 环境 |
|---|---|
| Python ≥3.8 + PyYAML、snakemake | base |
| fastp / SPAdes / BUSCO | `fastp` / `spades` / `busco` |
| GetOrganelle（含 animal_mt 库）/ MEANGS | `getorganelle` / `meangs` |
| jellyfish；R + genomescope、minpack.lm、argparse | PATH；`genome_size` |
| MitoM | `mitom`（`~/.local/bin`） |
| BLAST+、Biopython | PATH |

本机数据库：BUSCO `/media/deyuan/18T_A/database/busco/metazoa_odb12`；
18S-28S `/media/deyuan/18T_A/database/18-28s/`；nt `/media/deyuan/18T_A/database/blast_db/nt/nt`。

## 注意

- `mt-nr` 模式依赖 GetOrganelle 内置 animal_mt 种子，类群较远（如多毛类）时会报 "Too few seed reads"，改用 `all`（MEANGS 种子）或 `-s` 给近缘种。
- `fastp_spades_busco.py` 的 `--dryrun` 会返回非 0（把"无输出"记为失败），一键入口在 dryrun 时忽略。
- `fastp_spades_busco.py` 默认给 SPAdes 加 `--only-assembler`（跳过纠错），需要纠错时检查 `--spades_extra`。
- `genome_size.py` 默认生成 histo 后删除 `reads.jf`；`--keep_jf` 保留。低深度（<15×）时 GenomeScope 拟合不可靠，看 `model_fit` 和图。
