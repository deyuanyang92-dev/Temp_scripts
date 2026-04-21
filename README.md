# Batch GetOrganelle

`get.py` is a batch wrapper for GetOrganelle. It automates paired-end read
detection, multi-sample assembly, resume handling, result status reporting,
intermediate cleanup, and final `*.path_sequence.fasta` summary.

## Quick Start

```bash
python get.py -i reads/ -o out/ -F embplant_pt \
  --suffix_fq _R1.fq.gz _R2.fq.gz
```

This default mode runs:

```text
assembly -> clean -> summary
```

## Common Commands

Run assembly only:

```bash
python get.py assembly -i reads/ -o out/ -F embplant_pt
```

Clean intermediate files:

```bash
python get.py clean -o out/
```

Summarize final path sequences:

```bash
python get.py summary -o out/
```

Preview commands without running GetOrganelle:

```bash
python get.py assembly -i reads/ -o out/ -F embplant_pt --dry_run
```

## Notes

- Default completion requires non-empty `*.path_sequence.fasta`.
- Graph-only results are reported as `INCOMPLETE` unless `--accept_graph_only`
  is explicitly used.
- Failed or incomplete samples keep intermediate files by default for debugging.
- Resume safety is guarded by `batch_getorganelle_manifest.json`.

## Requirements

- Python 3
- GetOrganelle
- SPAdes
- Bowtie2
- BLAST+
- Biopython is required when converting GenBank files to FASTA.
