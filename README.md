# Temp_scripts

Small bioinformatics and utility scripts collected by topic.

## Directory Layout

```text
Temp_scripts/
  blast2metadata/
    blast2_metadatav2.12.py
    tests/
  getorganelle/
    batch_getorganelle.py
    README.md
  mitofinder-batch/
    mitofinder_pipeline_v2_7.py
    README.md
  Mitoz-annotate/
    batch_mitoz.py
    README.md
    docs/flowcharts/
    scripts/render_workflow_diagrams.py
```

## Available Scripts

### blast2metadata

Batch BLAST (local, `blastn -remote`, or Biopython qblast) -> TopN hits per
query -> GenBank metadata via Entrez -> merged CSV. Python or Snakemake
scheduling, `--resume`.

- TopN counts subjects (best HSP per query-subject), default sort `bitscore`.
- A failed BLAST batch fails the file (non-zero exit, no completion marker);
  `--resume` reruns only the failed batches.
- `--mode add-header` adds a header row to any existing BLAST tabular file
  (outfmt 6/7/10/20). Columns come from a `# Fields:` line or header row, from
  `--outfmt_blast "6 std staxids"`, or from the column count (12 = `std`), and
  are type-checked so a wrong guess is refused.
- `--header_style code|label|both` writes `pident`, `Percentage of identical
  matches`, or both as column headers (descriptions from `blastn -help`, 2.17).
- Behind an HTTP proxy `blastn -remote` failed or hung in testing; use
  `--online_engine biopython`.

Tests (real BLAST+): `mamba create -n blast -c conda-forge -c bioconda blast
biopython pandas psutil pytest snakemake-minimal`, then
`mamba run -n blast python -m pytest -q blast2metadata/tests`.

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
