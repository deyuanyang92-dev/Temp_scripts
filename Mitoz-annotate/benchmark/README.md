# mito_intron.py benchmark

Reproduce (paths are examples):

1. Positive set: NCBI nuccore `Metazoa[Organism] AND mitochondrion[Filter] AND complete genome[Title] AND intron[Feature key] AND refseq[Filter]` (85 records, 2026-10-09); negative set: 20 RefSeq annelid mitogenomes (Nereididae, Polynoidae, Sabellidae, Lumbricidae).
2. Split into `bench/fasta/<acc>.gb|.fasta` and `bench/records.tsv` (acc, set, genetic code, organism).
3. `python3 make_mitoz_like_input.py bench` -> `bench/sim/<acc>/` (introns removed, split genes cut to exon 1).
4. Ours: `mito_intron.py -i bench/sim/<acc> -o bench/ours_full/<acc> --genetic_code <code> --bin <env>/bin` (and `--exclude_own_species` -> `ours_loo`).
5. MFannot: `docker run --rm -v bench:/data -w /data/mfannot nbeck/mfannot mfannot -g <code> -o /data/mfannot/<acc>.mf /data/fasta/<acc>.fasta`.
6. `python3 compare.py bench` -> summary + `bench/compare.tsv`.

Truth = NCBI intron features + introns implied by joined CDS (gap >= 142 nt, not wrapping the origin) without an intron feature. A prediction matches a truth intron when it overlaps >= 50 % of it. Results: see ../README.md.

## eval_trna.py (mito_trna.py)

Records: `trna_sets.tsv` (85 intron-containing RefSeq mitogenomes, genetic code 4 as run; 40 annelid RefSeq
mitogenomes with 22 NCBI tRNAs, first 40 by accession from the Annelida `incl-refseq-mitogenome` dedup set).
Strip all features but `source`, run `mito_trna.py -i <acc_dir> -o out/<acc> --genetic_code <code>
--tools mitos2,arwen,aragorn,trnascan`, then `python3 eval_trna.py <ncbi_gb_dir> out <label>`.
