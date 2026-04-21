# MitoFinder Batch Pipeline

`mitofinder_pipeline_v2_7.py` 是一个面向 MitoFinder 的批量处理脚本，用于：

- 对 reads 批量组装；
- 对 reads 执行 MitoFinder 组装和注释；
- 对已有 contig/assembly fasta 执行 MitoFinder 注释；
- 汇总 MitoFinder 的 `*_Final_Results` 结果。

官方 MitoFinder 项目：<https://github.com/RemiAllio/MitoFinder>

## 1. 环境要求

### Python 依赖

脚本需要 Python 3 和 Biopython：

```bash
pip install biopython
```

### MitoFinder 可执行环境

正式运行 `assembly2summary` 或 `annotation2summary` 时，需要以下任一方式可用。

方式 1：`mitofinder` 命令在 `PATH` 中：

```bash
which mitofinder
```

方式 2：通过参数指定本地 MitoFinder 可执行脚本：

```bash
--mitofinder_script_path /path/to/mitofinder
```

方式 3：通过 Singularity 镜像运行：

```bash
--singularity_sif /path/to/mitofinder_v1.4.2.sif
```

官方 Singularity 镜像可按 MitoFinder README 的说明获取：

```bash
singularity pull --arch amd64 library://remiallio/default/mitofinder:v1.4.2
```

如果正式运行时找不到 MitoFinder，脚本会在任务池启动前退出，并提示以上参数，避免几十个样本同时报 `No such file or directory: 'mitofinder'`。

## 2. 子命令概览

查看帮助：

```bash
python3 mitofinder_pipeline_v2_7.py --help
```

主要子命令：

```text
only_assembly       仅批量组装 reads，不做 MitoFinder 注释
assembly2summary    reads -> MitoFinder 组装+注释 -> 汇总
annotation2summary  已有 contig fasta -> MitoFinder 注释 -> 汇总
snakemake           生成 Snakemake 工作流
analyze             分析输入目录结构，给出参数建议
```

## 3. annotation2summary：已有 fasta 注释

### 3.1 输入目录中有多个 fasta

如果目录中每个 `.fasta` 是一个样本：

```bash
python3 mitofinder_pipeline_v2_7.py annotation2summary \
  -i /data/contigs \
  -o ann \
  -r /path/to/reference.gb \
  --mitofinder_script_path /path/to/mitofinder
```

默认扫描：

```text
/data/contigs/*.fasta
```

### 3.2 单个样本 fasta

如果输入是单个样本的 fasta：

```bash
python3 mitofinder_pipeline_v2_7.py annotation2summary \
  -i sample.fasta \
  -o ann \
  -r /path/to/reference.gb \
  --mitofinder_script_path /path/to/mitofinder
```

脚本会生成一个 MitoFinder 注释任务：

```bash
mitofinder -j sample -a sample.fasta -r reference.gb ...
```

### 3.3 多样本合并 fasta 和重复 header

脚本支持 `rename_summary_all.fasta` 这类多样本合并 fasta。例如：

```text
>B2 topology=circular
...
>B2 topology=circular
...
>B17 topology=linear
...
>B17 topology=linear
...
```

处理规则：

- 使用 header 第一列作为样本名，例如 `B2`、`B17`；
- 同一样本名重复出现时，不去重、不覆盖、不丢弃；
- 认为重复记录是同一个样本的多条 contig；
- 每个样本拆成一个 fasta；
- 拆分后的重复 header 自动改名为 `sample_contigN`。

例如 `B2` 有 24 条记录，会拆成：

```text
ann/_split_input_fastas/rename_summary_all/B2.fasta
```

内部 header 为：

```text
>B2_contig1 topology=circular
>B2_contig2 topology=circular
...
>B2_contig24 topology=circular
```

然后 MitoFinder 会对 `B2` 这个样本运行：

```bash
mitofinder -j B2 -a ann/_split_input_fastas/rename_summary_all/B2.fasta -r reference.gb ...
```

这符合 MitoFinder 官方 `-a/--assembly` 的用法：一个 assembly fasta 可以包含一条或多条 contig。

运行示例：

```bash
python3 mitofinder_pipeline_v2_7.py annotation2summary \
  -i rename_summary_all.fasta \
  -o ann \
  -r /path/to/reference.gb \
  --mitofinder_script_path /path/to/mitofinder
```

如果需要强制按 header 第一列拆分单个 fasta：

```bash
--split_input_fasta_by_id
```

### 3.4 拆分映射表

拆分目录中会生成：

```text
split_manifest.tsv
```

字段：

```text
sample_name
new_contig_id
original_header
original_record_index
split_fasta
```

用途：

- 追溯原始 header；
- 确认重复 header 没有被覆盖；
- 检查每个样本拆出了多少条 contig。

例如统计 `B2` 条数：

```bash
awk -F '\t' '$1=="B2"{n++} END{print n+0}' \
  ann/_split_input_fastas/rename_summary_all/split_manifest.tsv
```

## 4. assembly2summary：reads 到注释汇总

目录中每个子目录是一个样本，脚本按 reads 后缀识别 R1/R2：

```bash
python3 mitofinder_pipeline_v2_7.py assembly2summary \
  --reads_dir /data/reads \
  -o ann \
  -r /path/to/reference.gb \
  --assembler megahit \
  --mitofinder_script_path /path/to/mitofinder
```

默认后缀：

```text
R1: *_1.clean.fq.gz
R2: *_2.clean.fq.gz
SE: *.fastq.gz
```

可通过参数修改：

```bash
--r1_suffix _R1.fastq.gz
--r2_suffix _R2.fastq.gz
--se_suffix .fastq.gz
```

## 5. only_assembly：仅组装

只运行外部组装器，不做 MitoFinder 注释：

```bash
python3 mitofinder_pipeline_v2_7.py only_assembly \
  --reads_dir /data/reads \
  -o assembly_out \
  --assembler megahit
```

支持：

```text
megahit
metaspades
idba
```

注意：当前脚本中 `idba` 自动转换交织 fasta 的流程未实现，建议优先使用 `megahit` 或 `metaspades`。

## 6. 常用参数

MitoFinder 核心参数：

```bash
-r, --genbank_ref FILE       参考线粒体 GenBank 文件
--genetic_code N             NCBI 遗传密码编号，默认 5
-p, --threads N              每个任务线程数，默认 4
-m, --memory GB              每个任务内存，默认 50
-t, --tool_param TOOL        tRNA 工具：mitfi/trnascan/arwen
```

并发参数：

```bash
--max_workers N              并行样本数
--task_timeout SEC           单样本超时时间
```

建议总线程数：

```text
总线程 = max_workers * threads
```

例如 16 核机器可设置：

```bash
--max_workers 4 --threads 4
```

## 7. 输出目录

如果指定：

```bash
-o ann
```

主要输出：

```text
ann/mitofinder_annotate/      annotation2summary 的 MitoFinder 结果
ann/mitofinder_assembly/      assembly2summary 的 MitoFinder 结果
ann/mitofinder_summary/       汇总结果
ann/_split_input_fastas/      合并 fasta 的自动拆分中间文件
```

汇总目录中包括：

```text
summary_mitofinder_results.tsv
gather_genes/
mt_contigs/
mitofinder_assembly_results/
```

## 8. Dry-run 检查

正式运行前建议先 dry-run：

```bash
python3 mitofinder_pipeline_v2_7.py annotation2summary \
  -i rename_summary_all.fasta \
  -o ann \
  -r /path/to/reference.gb \
  --dry_run
```

dry-run 会：

- 打印将要执行的 MitoFinder 命令；
- 检查样本拆分数量；
- 写出合并 fasta 拆分文件和 `split_manifest.tsv`；
- 不真正运行 MitoFinder。

## 9. 常见问题

### 找不到 mitofinder

错误：

```text
未找到 MitoFinder 可执行环境
```

解决方式：

```bash
which mitofinder
```

或：

```bash
--mitofinder_script_path /path/to/mitofinder
```

或：

```bash
--singularity_sif /path/to/mitofinder_v1.4.2.sif
```

### 重复 header 会不会丢序列？

不会。脚本不会按 header 去重，也不会用后面的记录覆盖前面的记录。同一样本重复记录会被保留为：

```text
sample_contig1
sample_contig2
sample_contig3
```

原始 header 可在 `split_manifest.tsv` 中追溯。

### 重复运行会不会覆盖拆分文件？

不会。若拆分目录已存在，脚本会创建：

```text
rename_summary_all_run2
rename_summary_all_run3
...
```

避免覆盖旧拆分结果。

