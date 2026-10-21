# Per-phylum 12S/16S profile HMMs (mito_rrna_hmm.py)

Built 2026-10-21 with `python3 mito_rrna_hmm.py build --gb refs/*.gb -o models/rrna_hmm`
(MAFFT --auto of <= 60 random reference rRNAs per phylum and gene, hmmbuild --dna; HMMER 3, MAFFT 7).

References: NCBI nuccore, per phylum `<phylum>[Organism] AND mitochondrion[Filter] AND (complete genome[Title]
OR complete sequence[Title]) AND 5000:60000[SLEN]`, up to 150 random records, with every accession of the
MitoM benchmark (benchmark/sets/mito_bench_2026-10.tsv, dev + test, RefSeq and source GenBank) removed
(knowledge/experiments/rrna_ref/fetch_refs.py); rRNAs 300-2500 nt with <= 1 % N. manifest.tsv lists the
phylum, gene, number of references, number aligned and median length. Phyla with < 3 references have no model.
