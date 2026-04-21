# Temp_scripts

Small bioinformatics and utility scripts collected by topic.

## Directory Layout

```text
Temp_scripts/
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
