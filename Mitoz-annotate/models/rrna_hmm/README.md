# Per-phylum 12S/16S profile HMMs (mito_rrna_hmm.py)

Rebuilt 2026-10-21 with `python3 mito_rrna_hmm.py build --gb refs/*.gb -o models/rrna_hmm --cpu 8`
(MAFFT --auto of <= 60 random reference rRNAs per phylum and gene (seed 1), hmmbuild --dna; HMMER 3, MAFFT 7).
41 models (21 phyla), 8 MB. Changes against the first build (same day): reference length cap raised from 2500
to 4000 nt (Porifera rrnL 80 -> 139 references, Cnidaria rrnL 129 -> 138 (incl. one relabelled), Placozoa rrnL 19 -> 22; before, 43 %
of the sponge 16S genes were dropped and the model came from the shorter ones); three references named
"16S-rRNA" were counted as 12S by the old name matcher and are now 16S (Rotifera 2, Cnidaria 1); manifest.tsv has
two more columns.

References: NCBI nuccore, per phylum `<phylum>[Organism] AND mitochondrion[Filter] AND (complete genome[Title]
OR complete sequence[Title]) AND 5000:60000[SLEN]`, up to 150 random records, with every accession of the
MitoM benchmark (benchmark/sets/mito_bench_2026-10.tsv, dev + test, RefSeq and source GenBank) removed
(fetched 2026-10-20 with knowledge/experiments/rrna_ref/fetch_refs.py; `scripts/fetch_rrna_refs.py` is its
reusable form and makes the same query, but was not used for this library). `refs_<phylum>.txt` lists the
accession.version of every record of that reference pool (21 files, one per phylum with a model); each model is
aligned from <= 60 of the pool's rRNAs. rRNAs of 300-4000 nt with <= 1 % N are used. Phyla with < 3 references
have no model (Entoprocta, Gastrotricha, Gnathostomulida, Kinorhyncha, Nematomorpha, Orthonectida, Phoronida,
Priapulida).

manifest.tsv columns: phylum, gene, n_ref (references passing the filters), n_used (aligned), median_len,
p5_len, p95_len (5th and 95th percentile of the reference lengths). mito_rrna_hmm.py accepts a hit when
it covers >= 40 % of the model (LENG in the .hmm header) and is <= 1.3 x p95_len; spans with both ends complete
must also be >= 0.7 x p5_len. `python3 mito_rrna_hmm.py build --gb refs/*.gb -o models/rrna_hmm --manifest-only`
recomputes manifest.tsv from the GenBank references without touching the .hmm files.

Known limits: Placozoa rrnL (22 features from 10 records) mixes the separate regions of a split 16S gene (p5 410,
p95 3008 nt, so its length window is very wide); models of phyla with 3-8 references (Onychophora, Brachiopoda,
Hemichordata, Chaetognatha, Tardigrada) are thin and only sensitivity, not specificity, was measured on them.
A model longer than the typical gene (Ctenophora rrnS: median 356 nt, p95 572) makes complete genes look
5'/3'-partial (hmm_from > 30) in the report.
