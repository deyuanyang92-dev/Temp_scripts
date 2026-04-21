# Batch GetOrganelle

`batch_getorganelle.py` is a batch wrapper for GetOrganelle. It automates paired-end read
detection, multi-sample assembly, resume handling, result status reporting,
intermediate cleanup, and final `*.path_sequence.fasta` summary.

## Requirements

- Python 3
- GetOrganelle
- SPAdes
- Bowtie2
- BLAST+
- Biopython is required when converting GenBank files to FASTA.

Before running real jobs, activate the environment where GetOrganelle is
installed and initialized.

## Quick Start

```bash
cd getorganelle
python batch_getorganelle.py -i reads/ -o out/ -F embplant_pt \
  --suffix_fq _R1.fq.gz _R2.fq.gz
```

This default mode runs:

```text
assembly -> clean -> summary
```

## Common Commands

Run assembly only:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_pt
```

Clean intermediate files:

```bash
python batch_getorganelle.py clean -o out/
```

Summarize final path sequences:

```bash
python batch_getorganelle.py summary -o out/
```

Preview commands without running GetOrganelle:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_pt --dry_run
```

Run only selected samples:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_pt \
  --sample_list samples.txt
```

Skip samples already confirmed elsewhere:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_pt \
  --finished_samples done.txt
```

Force rerun existing outputs, while still protecting `done.txt` samples:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_pt \
  --force --finished_samples done.txt
```

Use reads stored in selected subdirectories:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_pt \
  --subdir_list batch1 batch2 batch3
```

Use all reads for difficult or low-depth samples:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_pt --all_data
```

Plant mitochondrial assembly often needs more extension rounds:

```bash
python batch_getorganelle.py assembly -i reads/ -o out_mt/ -F embplant_mt -R 50
```

Custom target mode requires both seed and label references:

```bash
python batch_getorganelle.py assembly -i reads/ -o out_custom/ -F anonym \
  -s seed.fasta --genes label.fasta --max_extending_len 100 -P 0
```

Prepare seed/label files from GenBank or FASTA references:

```bash
python batch_getorganelle.py prep_db -o out/ --ref_files reference.gb --make_label
```

Use a GenBank reference directly during assembly:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_pt \
  --ref_gb reference.gb
```

## Input Naming

By default the script expects paired-end files like:

```text
SampleA_1.clean.fq.gz
SampleA_2.clean.fq.gz
```

For other naming patterns, set both suffixes:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_pt \
  --suffix_fq _R1.fq.gz _R2.fq.gz
```

## Output

Typical output:

```text
out/
  SampleA/
    get_org.log.txt
    *.path_sequence.fasta
    *.fastg
    batch_getorganelle_manifest.json
  summary_all.fasta
```

## Notes

- Default completion requires non-empty `*.path_sequence.fasta`.
- Graph-only results are reported as `INCOMPLETE` unless `--accept_graph_only`
  is explicitly used.
- Failed or incomplete samples keep intermediate files by default for debugging.
- Resume safety is guarded by `batch_getorganelle_manifest.json`.
- `--dry_run` previews commands without running GetOrganelle.
