# Changelog

## 3.0.0 (2026-10-22)

Packaged as `blast2metadata` (`pip install -e .` -> `blast2metadata` command, or
`python -m blast2metadata`, or `python run_blast2metadata.py` from the checkout).
The single 3300-line `blast2_metadatav2.12.py` is split into modules under
`src/blast2metadata/`; the command-line interface is unchanged except for the
additions below. Snakemake child jobs now run `python -m blast2metadata`.

New:
- `--engine diamond|mmseqs` for `-bm local`. DIAMOND (protein DB: blastp /
  blastx) and MMseqs2 (nucleotide or protein DB) write the same 14 columns as
  BLAST+, with `qcovs` computed from the union of HSPs per query-subject as
  BLAST does. `--blast_db` may be a ready `.dmnd` / MMseqs2 DB or a FASTA,
  which is converted once into `OUT/_engine_db/` (lock-protected, reused by
  `--resume` and parallel Snakemake jobs). Batch size defaults to 5000 and one
  job with all threads for these engines.
- `--engine diamond --blast_db nr`: DIAMOND >= 2.1.14 searches a BLAST v5 protein
  database directly (multi-volume `nr.pal` alias included), so the nr from
  `update_blastdb.pl --decompress --blastdb_version 5 nr` needs no makedb/.dmnd;
  bare names are looked up in `$BLASTDB` (DIAMOND itself ignores it), older
  DIAMOND or a nucleotide BLAST DB gives a clear error. MMseqs2 given a BLAST DB
  is pointed at `mmseqs databases NR|NT|core_nt`. `--doctor` reports whether the
  installed DIAMOND can read BLAST DBs.
- `blast2metadata-db` command (list / update / clean / env): lists the databases in
  `$BLASTDB`, compares them with NCBI's `<db>-nucl|prot-metadata.json` manifests,
  downloads a database (aria2c with md5 checksums and retry rounds, ascp, or
  update_blastdb.pl), extracts into a staging directory, verifies with
  `blastdbcmd -info` and swaps it in; cleans leftover tarballs. Replaces the
  ad-hoc aria2_nt.py / ascp_nt.py scripts. `--doctor` now also reports the
  downloaders and whether each `$BLASTDB` directory exists.
- `--query_gencode N` sets the genetic code of translated searches in one place
  (blastx `-query_gencode`, DIAMOND `--query-gencode`, MMseqs2 `--translation-table`).
- gzipped FASTA input (`*.fasta.gz`, output stem drops `.fasta.gz`).
- `--doctor` prints tool/package versions, proxy and BLASTDB; `--version`.
- `--topN`, `--batch_size`, `--max_target_seqs`, `--retries` must be >= 1
  (rejected at parse time).

Robustness:
- A failed Entrez efetch batch (typically one invalid accession -> HTTP 400 for
  the whole batch) is split in halves and retried, so the valid accessions still
  get metadata.
- Ctrl-C cancels the batches not yet started instead of waiting for all of them.
- Biopython qblast now sends the `--email` / `--entrez_tool` (previously only
  Entrez did).
- Duplicate sequence IDs in a query FASTA are reported (their hits pool in TopN).
- Engine executables/DB paths are made absolute before Snakemake children start.

Tests: real-tool tests for DIAMOND and MMseqs2 (skipped when not installed; set
`DIAMOND_EXE` / `MMSEQS_EXE`), gzip input, Entrez batch splitting, entry points.
The outfmt-20 header test is skipped on BLAST+ < 2.13.

## 2.12 and earlier

See the git history of `blast2metadata/blast2_metadatav2.12.py`; the patch
notes that used to head that file (Fix-1 ... Fix-J, Codex 修改, 2026-10-10/11
changes) describe the TopN-by-subject fix, `.part` outputs for safe `--resume`,
the split Entrez/BLAST rate limiters, `--args_json` Snakemake children,
`--mode add-header` and the BLAST field dictionary.
