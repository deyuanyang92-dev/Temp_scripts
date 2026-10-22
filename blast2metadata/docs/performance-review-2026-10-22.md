# BLAST 与后处理评估（2026-10-22）

结论：原版不能称为“性能最优、最稳健”。其默认 bitscore 排序、按 subject 去除重复 HSP、计算 HSP 并集覆盖度的方向合理，但实现中存在资源超额、内存放大、缓存误用和结果语义混淆。本次已修复可复现的问题，并保留以下性能与科学解释边界。

## 实测性能

环境：Linux 6.8，Python 3.12.9、Biopython 1.85，本机 `blastn 2.12.0+`。同一机器上独立进程重复 3 次；报告中位数，峰值内存为整个 Python 进程的 Linux `ru_maxrss`，包含导入的库。共享机器负载会影响时间，不能把单次最好成绩当作保证。

| 测量 | 原版/原配置 | 优化后 | 结果 |
|---|---:|---:|---|
| 30 万行解析 + Top10 | 10.851 秒 | 4.693 秒 | 约 2.31 倍速度 |
| 同一任务峰值 RSS | 368.6 MiB | 76.5 MiB | 降低约 79% |
| 本地小库 1500 条查询 | 50 条/批：0.688 秒 | 500 条/批：0.246 秒 | 约 2.80 倍速度 |

TopN 基准：300 个 query、每 query 100 个 subject、每对 10 个 HSP，共 300000 行；按 bitscore 取 Top10，输出 3000 行。新旧输出按键序列化后逐字一致。基准刻意避免本次修正的版本 ID 合并和同分排序差异。数据是合成表，不代表所有 query/subject 分布。

BLAST 基准：使用仓库 `tests/data/mito_db.fasta`，把 6 条约 600 bp 查询复制成 1500 个不同 ID，线程预算不超过 4。50 条/批使用 4 个进程、每进程 1 线程；500 条/批实际只有 3 批，因此使用 3 个进程、每进程 1 线程。计时不含建库和准备 FASTA，所有方案的 HSP 行集合完全相同。500 条/批配单进程的中位数为 0.677 秒，1500 条/批配单进程为 0.513 秒：**批次越大、单进程线程越多，不一定越快。**

重现命令：

```bash
python benchmarks/benchmark_topn.py
# old_blasttab.py 是优化前 src/blast2metadata/blasttab.py 的副本：
python benchmarks/benchmark_topn.py --baseline-module /path/to/old_blasttab.py
python benchmarks/benchmark_blast_batches.py
```

原始测量见 `benchmark-topn.jsonl` 和 `benchmark-blast.jsonl`。没有对完整 nt/nr 进行生产规模搜索，也没有改变 word size、task、E-value 或屏蔽策略来换取速度。

## BLAST 调度与稳健性

- **线程预算**：原版 `ceil(total/jobs)` 在 10 线程、4 任务时实际分配 12 线程；手动指定每任务线程也不受上限约束。现改为向下分配，必要时减少并行任务，保证每次运行的分配乘积不超预算。
- **批次**：本地 BLAST 默认 50 → 500（2026-10-22）→ **5000**（2026-10-23，依据见下节真实 nt 测量）；在线仍为 50，DIAMOND/MMseqs2 为 5000。用户显式 `--batch_size` 保持优先。
- **续跑**：输入内容 SHA-256、搜索参数、程序文件状态和本地数据库文件状态构成搜索签名；后处理另有参数签名。修改 TopN、过滤阈值或表头时复用搜索、重算结果；改变输入或搜索参数会重跑。零命中成功批次也有缓存凭据。
- **原子输出**：TSV/CSV/完成标记写完 `.part` 后替换；重试清理上一轮半成品。qblast 格式化中途失败时撤销本条 query 已写入的 HSP，避免重试后重复。
- **失败处理**：支持 `--blast_timeout`，默认 0 保留不限时行为；有缺失 metadata 时保存部分结果和缺失 accession 清单，非零退出、不写完成标记，允许续跑补齐。
- **输入与参数**：重复 query ID、空序列提前失败；不能通过 extra 参数覆盖线程和受管理的输入/输出参数。翻译搜索按 subject 类型选择默认 nt/nr，显式程序绝对路径也能识别。
- **Snakemake**：参数 JSON 纳入任务输入，防止参数改变后调度器仅因 CSV 已存在而跳过签名检查。

### 2026-10-23 补充：真实 nt 测量

小库基准不能代表 nt。改用本机 nt 前 10 卷（28 GB，已预热进页缓存），blastn 2.12.0（默认 megablast），总线程 8，600 bp COI 突变副本（3% 替换，每条不同）。测量时本机负载约 55/72（其他任务），只测 1 轮，绝对值有噪声，但差距远大于噪声。所有配置 HSP 集合完全一致。原始数据 `benchmark-nt-passes.jsonl`，脚本 `benchmarks/benchmark_nt_passes.py`。

| 配置 | 秒 |
|---|---:|
| 6 条，1 批，1 进程 × 8 线程 | 13.7 |
| 96 条，1 批，1 进程 × 8 线程 | 20.0 |
| 960 条，1 批，1 进程 × 8 线程 | **66.9** |
| 960 条，2 批，2 进程 × 4 线程 | 77.8 |
| 960 条，20 批（≈旧默认 50/批），4 进程 × 2 线程 | 211.9 |
| 960 条，1 批，1 进程 × 8 线程，`-mt_mode 1` | 407.6 |

结论：

- 每批约 13 秒是**扫库的固定开销**，与 query 数几乎无关（6 条与 96 条只差 6 秒）。批越多，扫库次数越多；20 小批比 1 大批慢 3.2 倍。完整 nt 的 `.nsq` 共 800 GB，超过本机 503 GB 内存，每多扫一遍还要多读一遍磁盘，差距只会更大。因此本地 BLAST 默认批大小改为 5000。megablast 内部一次最多拼接约 5 Mbp query，5000 × 600 bp 仍在一遍之内；更长的 query 由 BLAST 自己内部分遍。代价：单批失败时 `--resume` 重跑的序列更多。
- `-mt_mode 1`（按 query 分线程）在 nt 上慢 6 倍：每个线程各扫一遍库。**不采用**，保持默认 mt_mode 0（按库分线程）。
- 2 进程 × 4 线程比 1 进程 × 8 线程慢 16%。但本机负载高，未能测 16–72 线程下单进程 BLAST 的扩展性，所以**未改** `--parallel_jobs` 默认（≤4）。批大小 5000 下，5000 条以内的输入本来就只有 1 批，自动成为 1 进程占满全部线程。
- Snakemake 多样本时每个样本至少扫一遍库；样本多而每个样本序列少时，把样本合并成一个 FASTA 跑一次，再按 query ID 拆结果，会比并行跑多个样本快。本次未改 Snakemake 调度。

NCBI 从 BLAST+ 2.15.0 起提供自动线程模型选择，适合保留程序自动决策；本机 2.12.0 不具备该版本行为。查询文件大小、数据库大小和搜索 task 都影响最合适的线程模型。[NCBI 线程模型文档](https://www.ncbi.nlm.nih.gov/books/NBK571452/)

### 2026-10-23 补充：在线 BLAST（qblast）

原实现每条序列单独 `NCBIWWW.qblast`。Biopython 1.85 提交后先等 20 s 再取结果，之后每 60 s 一次（NCBI 规定同一 RID 每分钟最多查询一次），再加每条固定 sleep 12 s，所以耗时按“条数 × 至少一次等待”计算。NCBI URL API 支持一次提交多条 query，XML 中每条一个 Iteration。

实测（真实 NCBI nt，经 HTTP 代理，数据 `benchmark-online-qblast.jsonl`）：

| 方式 | 条数 | 秒 | 秒/条 |
|---|---:|---:|---:|
| 旧：逐条提交 + 12 s | 3 | 246.3 | 82.1 |
| 新：一批一次，URL API 默认（普通 blastn） | 10 | 125.2 | 12.5 |
| 新：一批一次，MEGABLAST=on | 10 | 27.7 | 2.8 |
| 新：完整 CLI blast2topN，一批 50 条 | 50 | 72.1 | 1.4 |

- 逐条与成批的命中对完全相同。
- URL API 不加 MEGABLAST=on 时是普通 blastn（2/−3 计分），与 `-remote`、本地 BLAST+、NCBI 网页默认的 megablast（1/−2）不同：最佳命中 pident 一致，但 bitscore 不同，同分较多的 query TopN 成员会变（10 条中 2 条前 10 重合仅 2、4 个）。现默认 megablast，使三条路径可比。
- 提交间隔按 NCBI 规定改为 ≥10 s；原 3.1 s 限速只因每次另有 12 s 固定等待才未超速。
- 单个样本只测一次，NCBI 排队时间随时段变化；倍数可信，绝对秒数不保证。
- 未改 `blastn -remote`：它本来就按批提交，轮询由 BLAST+ 控制；只是去掉了默认 12 s 额外等待。本机代理下 `-remote` 不可用，未实测。

## 后处理是否合理

保留默认 **bitscore → evalue → pident**：相比只按 identity，能避免很短的 100% 命中轻易压过更长的高相似命中。每个 subject 只保留当前排序下最好的 HSP，避免一个 subject 的多段比对占满 TopN。这里的分数是代表 HSP 分数，**不是 subject 全部 HSP 的综合证据评分**。

修正事项：

1. 按完整 `sseqid` 去重。原来截断第一个点后的所有内容，会合并不同 accession 版本和 `contig.partA/partB`。同分增加 subject ID/坐标/title 次序，合并批次按批次名排序，保证常见同分情况下结果不受 worker 完成顺序影响。
2. TopN 流式消费任意顺序的 HSP，只保存每 query 的 N 个候选。候选存储 O(query 数 × N)，无需全量表和全量排序。极大的 TopN 仍会增加时间和内存；当前有序列表插入为 O(N)。
3. 正确读取标准 12 列、14 列、outfmt 7、CSV，以及本软件写出的带引号/换行 TSV。缺列、非法数字不再当作“没有命中”。缺失 `qcovs` 时，要求覆盖度的排序或过滤明确失败。
4. 只开启 identity 阈值时，不再因为 coverage 列缺失而错误过滤该行。自定义排序限制为数值列，避免按 `stitle` 等文本列排序时所有值都被视作缺失。
5. accession 提取使用集合去重，避免大表逐项扫描列表的二次复杂度；本地 ID 中非版本含义的点不会被删除。

`qcovs` 是同一 subject 所有 HSP 的 query 覆盖度；`pident` 是单条 HSP 的一致率，两者粒度不同。[NCBI 字段定义](https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.options_common_to_all_blast/)

因此 `pident × qcovs / 100` 只作为排序启发式，不是精确“整条 query 一致率”：多 HSP 的 identity 不同、重叠区域、gap 和覆盖度取整都会影响解释。阈值作用于每个 subject 的**代表 HSP**（2026-10-23 修正，原先先逐 HSP 过滤再选代表：主 HSP 900 bp/85% 未达 97% 的 subject，可借一段 40 bp/100% 的短 HSP 通过 `--min_pident 97`，且 subject 级 `qcovs` 包含那段 85% 比对，也能通过 `--min_qcovs`，在库中无同种序列时报出假的高一致率命中；现与 NCBI 网页一致，Per. Ident 取 subject 最佳 HSP）。`qcovs` 沿用搜索报告值，不会按过滤后的 HSP 重算。没有直接把所有 HSP 的 bitscore 相加，因为重叠区域可能被重复计分。

TopN 表仍只是候选注释，不能自动等同物种鉴定结论。搜索阶段的 `max_target_seqs` 限制候选集合，事后换排序无法恢复未返回的 subject；低覆盖、多个同分物种、数据库注释错误、嵌合或污染需要结合原始比对和领域标准判断。

## 验证与仍存在的边界

回归覆盖真实 BLAST 本地搜索、失败批次续跑、输入/参数变化、零命中、TopN 与独立全排序算法对照、标准表格、坏行、特殊 title 往返、超时、metadata 缺失重试、翻译搜索库选择和 Snakemake 调度。最终 `python -m pytest -q tests`：**55 passed, 5 skipped**（55.37 秒）；`ruff check` 与 `git diff --check` 均通过。5 项跳过涉及本机未满足的 BLAST 字段版本、可选引擎或联网条件；未实际调用公网 BLAST/Entrez 做压力测试。Snakemake 测试需要用户缓存目录写权限，已在获准后完成。

需知的限制：

- `blast-only` / `blast2meta` 仍会加载完整结果；本次主要降低使用 TopN 路径的内存。DIAMOND/MMseqs2 的覆盖度重写仍在批次内缓存行。
- 大型数据库缓存验证使用路径、大小、mtime，避免散列 TB 级数据；不能发现保留大小和时间戳的原地修改，也不能完整追踪任意 BLAST alias 指向的外部数据库。在线数据库变化无法自动追踪。需要强制刷新时去掉 `--resume` 或使用新输出目录。
- 原有 DIAMOND/MMseqs2 从 FASTA 自动建库缓存仍按文件 stem 命名；同名参考替换/不同目录同名参考应使用新的输出目录，不能依赖旧 `_engine_db`。本次没有声称已完成这些引擎的全量缓存重构。
- 不能让多个独立任务同时写同一输出 stem；签名机制不是跨进程写锁。
- `--blast_timeout` 限制外部搜索直接子进程，不涵盖 qblast/Entrez 网络内部等待，也不承诺清理任意包装脚本启动的孙进程。Ctrl-C 仍可能等待已运行的 worker，建议长任务显式设置超时。
- metadata 查询仍按基础 accession 归并，可能获取当前记录版本；需要严格复现历史 accession.version 注释时，不能将其当作历史版本归档。
- 默认 `blastn` 程序通常运行 megablast task；是否改用 `-task blastn`、短序列 task 或不同遗传密码应按生物学问题选择，不能统一为提速而降低敏感性。

兼容性变化：旧版无签名缓存首次续跑会重算；重复 FASTA ID 和损坏表格从警告/静默跳过变为失败；metadata 未抓全从“成功”改为可续跑的失败状态；不同 subject 版本不再合并。关键改动旁已添加中文原因说明。
