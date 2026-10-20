# 细胞器注释知识库 / Organelle annotation knowledge base

用途：各注释软件的算法、源码关键位置、默认参数、已知 bug、类群局限，带出处，供优化 Mitoz-annotate 流程时查证，避免重复学习。
规则：每条结论注明出处（`path:line`、DOI、URL）、软件版本与查阅日期；没有依据的标 `未确证`。检查：`python3 benchmark/check_knowledge.py knowledge/`。

## tools/

| 文件 | 范围 |
|---|---|
| tools/mitoz.md | MitoZ 3.6 annotate: PCG (GeneWise/HMM), rRNA (Infernal global), tRNA (MiTFi), stop/start revision, clades |
| tools/mitfi.md | MiTFi tRNA finder inside MitoZ/MITOS2: options, naming, single-core |
| tools/mitos2.md | MITOS2 runmitos: BLAST/HMM proteins, Infernal ncRNA, naming (trnM_0, trnQ(---)), sensitive pass, refseq89m (RefSeq 89 = 2018-07-09) |
| tools/arwen.md | ARWEN metazoan mt tRNA finder (same engine as ARAGORN -mt) |
| tools/aragorn.md | ARAGORN tRNA/tmRNA finder; -mt / -mtx flags |
| tools/trnascan-se.md | tRNAscan-SE 2.x organellar mode -O, -M vert, genetic-code flag -g |
| tools/mfannot.md | MFannot intron workflow (exonerate, RNAweasel), .mfannotrc defaults; basis of mito_intron.py |
| tools/mitofinder.md | MitoFinder: reference-based annotation, py2.7 bioconda, needs related reference GenBank |
| tools/table2asn.md | NCBI table2asn validator: messages relevant to mitogenomes, transl_except for T/TA stops; experiments in experiments/table2asn |

## topics/

| 文件 | 范围 |
|---|---|

## 未确证问题 / Open questions

- refseq89m 实际基于哪个 RefSeq 版本（文件名指 89；Zenodo 描述写 63 与 81）— mitos2.md
- MitoZ 线程数依赖的缺基因原因；-onlycutoff 是否真导致 tRNA 漏注 — mitoz.md / mitfi.md
- NCBI 是否对动物线粒体 intron 检查剪接共有序列 — table2asn.md
- MitoFinder 对不完整终止 T/TA 的处理 — mitofinder.md
- ARWEN/ARAGORN 在无臂线虫 tRNA 上的表现 — arwen.md

## 已发现的本流程问题（待 Phase 2 修）/ Pipeline issues found

- mito_trna.py:133 tRNAscan-SE `-O` 未加 `-g`（遗传密码）→ Trp(TCA) 被判为 Sup 并被丢弃 — trnascan-se.md
- mito_trna.py:160 MITOS2 tRNA-only 未加 `--orih 0` → 多余 OH BLASTN（不影响 tRNA，费时）— mitos2.md
- mito_trna.py ARWEN 用 `-mtx`、ARAGORN 用 `-mt`：同一引擎不同参数，ARAGORN FP 更多 — arwen.md / aragorn.md
- ARWEN/ARAGORN 反密码子 `???`（6/8 nt 环）被解析丢弃 — arwen.md
- mito_intron.py MAX_INTRON 5000 vs MFannot 默认 maxintronsize 3000 — mfannot.md（需判断哪个对）

## benchmark_findings.md
每轮 benchmark 暴露的问题 → 原因 → 修复 commit。
