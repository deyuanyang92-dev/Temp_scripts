# Temp_scripts

Small bioinformatics and utility scripts collected by topic.

## Directory Layout

```text
Temp_scripts/
  getorganelle/
    batch_getorganelle.py
    README.md
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
- safer status reporting for `OK`, `SKIP`, `INCOMPLETE`, and `FAILED`.

See [getorganelle/README.md](getorganelle/README.md) for usage.
