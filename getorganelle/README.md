# Batch GetOrganelle

`batch_getorganelle.py` is a research-oriented batch wrapper for
GetOrganelle. It is designed for daily organelle assembly work, especially
batch assembly of plant mitochondrial genomes (`embplant_mt`) and nuclear
ribosomal regions (`embplant_nr` / `fungus_nr`).

The goal is not to replace GetOrganelle. The goal is to make the GetOrganelle
workflow easier to repeat across many samples: find paired reads, prepare
seed/label references, run samples in parallel, resume safely, iterate when
needed, clean successful runs, and summarize final `*.path_sequence.fasta`
files.

Official references:

- GetOrganelle: <https://github.com/Kinggerm/GetOrganelle>
- Usage: <https://github.com/Kinggerm/GetOrganelle/wiki/Usage>
- FAQ: <https://github.com/Kinggerm/GetOrganelle/wiki/FAQ>

## Requirements

- Python 3
- GetOrganelle
- SPAdes
- Bowtie2
- BLAST+
- Biopython for GenBank-to-FASTA conversion

Before running real jobs, activate the environment where GetOrganelle is
installed and its database has been initialized.

## Main Features

### 1. Batch GetOrganelle assembly

The script wraps `get_organelle_from_reads.py` for many paired-end samples.
It detects read pairs by suffix, builds one output directory per sample, runs
samples in parallel, and reports `OK`, `SKIP`, `INCOMPLETE`, and `FAILED`.

### 2. Batch mitochondrial genome assembly

The primary use case is batch mitochondrial assembly:

```bash
python batch_getorganelle.py assembly -i reads/ -o mt_out/ \
  -F embplant_mt -R 50 \
  --suffix_fq _R1.fq.gz _R2.fq.gz
```

Plant mitogenomes can be structurally complex. The script keeps graph-only or
broken results as `INCOMPLETE` by default, so they are not silently mixed into
the final sequence summary.

### 3. `-s` seed support: FASTA or GenBank

You can provide a custom seed with `-s/--seed`.

Supported seed formats:

- `.fasta`
- `.fa`
- `.fna`
- `.fas`
- `.gb`
- `.gbk`
- `.gbf`
- `.gbff`
- `.genbank`

When a GenBank-format file is provided to `-s`, the script automatically
converts it to FASTA with Biopython before passing it to GetOrganelle:

```bash
python batch_getorganelle.py assembly -i reads/ -o mt_out/ \
  -F embplant_mt -R 50 \
  -s reference_mt.gb
```

This produces a converted seed under:

```text
mt_out/db/
```

### 4. `--ref_gb`: automatically prepare seed and `--genes`

For annotated GenBank references, use `--ref_gb`. The script can:

- convert the GenBank sequence into seed FASTA;
- call GetOrganelle's `get_annotated_regions_from_gb.py`;
- extract annotation regions as a label FASTA;
- pass that label FASTA to GetOrganelle as `--genes`.

Example for mitochondrial assembly:

```bash
python batch_getorganelle.py assembly -i reads/ -o mt_out/ \
  -F embplant_mt -R 50 \
  --ref_gb reference_mt.gb
```

Label extraction type defaults to `CDS`. You can change it:

```bash
python batch_getorganelle.py assembly -i reads/ -o mt_out/ \
  -F embplant_mt -R 50 \
  --ref_gb reference_mt.gb \
  --ref_gb_region gene
```

Important: according to the GetOrganelle FAQ, using seed/label databases does
not make the assembly reference-guided. The seed helps recruit reads, and the
label helps graph filtering/disentangling. The assembly itself remains de novo.

### 5. Batch nr assembly with `--ref_gb`

The script also supports batch assembly of nuclear ribosomal regions:

```bash
python batch_getorganelle.py assembly -i reads/ -o nr_out/ \
  -F embplant_nr -R 10 -k 35,85,105,115 \
  --ref_gb reference_nr.gb \
  --ref_gb_region rRNA
```

For fungal nr:

```bash
python batch_getorganelle.py assembly -i reads/ -o nr_out/ \
  -F fungus_nr -R 10 -k 35,85,105,115 \
  --ref_gb reference_nr.gb \
  --ref_gb_region rRNA
```

nr results may be circular, linear, or multiple sequences depending on repeat
structure and concerted evolution. Inspect `*.fastg`/`*.gfa` and logs when
results are marked `INCOMPLETE`.

### 6. Iterative mt/nr assembly

For difficult mt or nr datasets, use `--iter_rounds`. Each completed round
uses the previous round's `*.path_sequence.fasta` as the seed for the next
round:

```bash
python batch_getorganelle.py assembly -i reads/ -o mt_iter_out/ \
  -F embplant_mt -R 50 \
  --ref_gb reference_mt.gb \
  --iter_rounds 3
```

For nr:

```bash
python batch_getorganelle.py assembly -i reads/ -o nr_iter_out/ \
  -F embplant_nr -R 10 -k 35,85,105,115 \
  --ref_gb reference_nr.gb \
  --ref_gb_region rRNA \
  --iter_rounds 3
```

In iterative mode, summary only uses the highest completed `iter_N` directory
for each sample, avoiding duplicate reporting of earlier rounds.

## Quick Start

```bash
cd getorganelle
python batch_getorganelle.py -i reads/ -o out/ -F embplant_mt -R 50 \
  --suffix_fq _R1.fq.gz _R2.fq.gz
```

Default mode runs:

```text
assembly -> clean -> summary
```

For safer testing, preview commands first:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_mt -R 50 \
  --dry_run
```

## Input Layout

By default, the script expects paired-end files like:

```text
SampleA_1.clean.fq.gz
SampleA_2.clean.fq.gz
SampleB_1.clean.fq.gz
SampleB_2.clean.fq.gz
```

For other file names, set suffixes explicitly:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_mt \
  --suffix_fq _R1.fq.gz _R2.fq.gz
```

If reads are grouped under subdirectories:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_mt \
  --subdir_list batch1 batch2 batch3
```

## Sample Selection

Run only selected samples:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_mt \
  --sample_list samples.txt
```

`samples.txt` can contain sample names:

```text
SampleA
SampleB
```

It can also contain absolute R1 paths. R2 is inferred from `--suffix_fq`:

```text
/data/project/SampleC_R1.fq.gz
```

Protect samples that should never be rerun:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_mt \
  --force --finished_samples done.txt
```

## Output

Typical output:

```text
out/
  SampleA/
    get_org.log.txt
    *.path_sequence.fasta
    *.fastg
    *.selected_graph.gfa
    batch_getorganelle_manifest.json
  summary_all.fasta
```

For iterative runs:

```text
out/
  SampleA/
    iter_0/
    iter_1/
    iter_2/
    iter_0_seed.fasta
    iter_1_seed.fasta
  summary_all.fasta
```

## Result Status

- `OK`: final `*.path_sequence.fasta` exists.
- `SKIP`: sample was already completed or listed in `--finished_samples`.
- `INCOMPLETE`: GetOrganelle returned successfully but no final path sequence
  was found. Inspect graph and log files.
- `FAILED`: command failed or resume safety checks failed.

By default, only `OK` samples are cleaned. `INCOMPLETE` and `FAILED` outputs are
kept for debugging.

## Useful Commands

Clean intermediate files after assembly:

```bash
python batch_getorganelle.py clean -o out/
```

Summarize final path sequences:

```bash
python batch_getorganelle.py summary -o out/
```

Prepare reusable seed/label files:

```bash
python batch_getorganelle.py prep_db -o out/ \
  --ref_files reference.gb reference2.gbf \
  --make_label
```

Use all reads for difficult or low-depth samples:

```bash
python batch_getorganelle.py assembly -i reads/ -o out/ -F embplant_mt \
  --all_data
```

## Notes

- `--ref_gb` has lower priority than explicitly provided `-s` and `--genes`.
- `-F anonym` requires both seed and label references.
- `--continue` is protected by `batch_getorganelle_manifest.json`; if reads or
  core parameters changed, the script refuses unsafe resume unless
  `--ignore_manifest_mismatch` is used.
- Graph-only results are not considered final by default. Use
  `--accept_graph_only` only after manual inspection.
