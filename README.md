# Temp_scripts

Small bioinformatics and utility scripts collected by topic.

## Directory Layout

```text
Temp_scripts/
  blast2metadata/
    pyproject.toml
    run_blast2metadata.py
    src/blast2metadata/
    tests/
    README.md
  getorganelle/
    batch_getorganelle.py
    README.md
  mitofinder-batch/
    mitofinder_pipeline_v2_7.py
    README.md
  wgs_skim_pipeline/
    run_pipeline.py
    config.example.yaml
    scripts/
    README.md
  Mitoz-annotate/
    batch_mitoz.py
    README.md
    docs/flowcharts/
    scripts/render_workflow_diagrams.py
```

## Available Scripts

### blast2metadata

Installable package (`pip install -e blast2metadata` -> `blast2metadata`
command; or `python blast2metadata/run_blast2metadata.py`). Batch search
(local BLAST+, `blastn -remote`, Biopython qblast, DIAMOND or MMseqs2) -> TopN
hits per query -> GenBank metadata via Entrez -> merged CSV. Python or
Snakemake scheduling, `--resume`, gzipped FASTA input.

- `--engine diamond|mmseqs` (with `-bm local`) uses DIAMOND (protein DB only:
  blastp/blastx) or MMseqs2 (nucleotide or protein DB) and produces the same 14
  BLAST columns. DIAMOND >= 2.1.14 reads the NCBI BLAST nr directly
  (`--blast_db nr`, found via `$BLASTDB`), no makedb needed; MMseqs2 uses
  `mmseqs databases NR|NT`; both also accept a FASTA, converted once to
  `OUT/_engine_db/`.
- TopN counts subjects (best HSP per query-subject), default sort `bitscore`.
  `--sort_by` takes a preset (`bitscore`, `identity`, `coverage`, `evalue`,
  `pident_x_qcovs` = pident*qcovs/100) or columns in order, e.g. `qcovs,pident`;
  `--min_pident` / `--min_qcovs` filter before TopN. For species ID e.g.
  `--min_qcovs 80 --sort_by pident_x_qcovs`.
- Pick best hits by ranking afterwards, not with `-max_target_seqs 1`: that
  limit is applied early in the search and need not return the best hit
  (Shah et al. 2018/2019, Bioinformatics; NCBI reply by Madden et al.).
- A failed BLAST batch fails the file (non-zero exit, no completion marker);
  `--resume` reruns only the failed batches. A failed Entrez batch is split and
  retried so one bad accession does not lose the others' metadata.
- `--mode add-header` adds a header row to any existing BLAST tabular file
  (outfmt 6/7/10/20); `--header_style code|label|both`.
- `blast2metadata-db list --remote | update nt -o DIR | clean | env` manages the
  local NCBI databases (manifest comparison, md5-verified aria2/ascp/
  update_blastdb.pl download, staged extraction and swap).
- `--doctor` reports tools and packages found; behind an HTTP proxy
  `blastn -remote` fails or hangs, use `--online_engine biopython`.

Tests (real BLAST+, optional DIAMOND/MMseqs2): `python -m pytest -q
blast2metadata/tests` (set `DIAMOND_EXE` / `MMSEQS_EXE` to include the engine
tests). See [blast2metadata/README.md](blast2metadata/README.md).

### wgs_skim_pipeline

WGS 浅测序流程：fastp（WGS/RNA）→ MEANGS 种子 + GetOrganelle 线粒体 / 18S-28S → SPAdes（多种
组装模式）→ BUSCO → jellyfish + GenomeScope2 基因组大小 → MitoM 线粒体注释 → BLAST 物种验证
（在线或本地 nt）。每个工具可单独运行；`run_pipeline.py -c config.yaml` 一键串联，断点续跑。

See [wgs_skim_pipeline/README.md](wgs_skim_pipeline/README.md) for usage.

### getorganelle

Research-oriented batch wrapper for GetOrganelle. The main purpose is to make
daily mitochondrial (`embplant_mt`) and nuclear ribosomal (`embplant_nr` /
`fungus_nr`) assembly work easier across many paired-end samples.

Main capabilities:

- batch assembly with `get_organelle_from_reads.py`;
- custom seed support via `-s` using FASTA or GenBank files;
- automatic GenBank-to-FASTA conversion for `.gb`, `.gbk`, `.gbf`, `.gbff`,
  and `.genbank`;
- `--ref_gb` support to prepare seed and `--genes` label files using
  GetOrganelle utilities;
- iterative mt/nr assembly with previous round results used as the next seed;
- clean mode for GetOrganelle intermediate files;
- summary mode for final `*.path_sequence.fasta` files;
- safer status reporting for `OK`, `SKIP`, `INCOMPLETE`, and `FAILED`.

See [getorganelle/README.md](getorganelle/README.md) for usage.

### mitofinder-batch

Batch MitoFinder pipeline for reads-to-annotation and existing-contig
annotation workflows. It supports MitoFinder execution, summary generation,
split handling for merged FASTA inputs, and batch-oriented output reporting.

See [mitofinder-batch/README.md](mitofinder-batch/README.md) for usage.

### Mitoz-annotate

Batch MitoZ annotation and GenBank post-processing pipeline. It supports FASTA
or GenBank input, internal ID mapping, optional two-pass annotation,
gene-based reorientation, original metadata transfer, `/organism=` restoration,
thread or Snakemake scheduling, bilingual documentation, and UTF-8 flowcharts
with CJK font support.

See [Mitoz-annotate/README.md](Mitoz-annotate/README.md) for usage.
