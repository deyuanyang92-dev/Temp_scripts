# ARAGORN

Every bullet outside "What it does" / "Relevance to our pipeline" ends with a source:
`(path/in/source.py:123)`, `(doi:10.xxxx/...)`, `(https://...)`, or is marked 未确证.
Check: `python3 benchmark/check_knowledge.py knowledge/`

Path convention: `aragorn1.2.41.c:N` = line N of the source file aragorn1.2.41.c (sha256 92a31cc5c0b0ad16d4d7b01991989b775f07d2815df135fe6e3eab88f5e97f4a, downloaded 2026-10-20 from the pyaragorn repository); `arwen1.2.3.c:N` = line N of the ARWEN source (see arwen.md).

## Version checked
- ARAGORN v1.2.41 (source header "Version 1.2.41 September 2nd, 2020"), conda package aragorn-1.2.41-h7b50bb2_5 from bioconda in env `mitotrna`, binary /root/miniconda3/envs/mitotrna/bin/aragorn (478,176 bytes); `aragorn -h` prints "ARAGORN v1.2.41", checked 2026-10-20 (aragorn1.2.41.c:30) (https://conda.anaconda.org/bioconda/linux-64/aragorn-1.2.41-h7b50bb2_5.tar.bz2)
- The source file has the same sha256 as the one in the bioconda recipe (92a31cc5...f97f4a), so the line numbers below describe the compiled program; the recipe's source URL is www.ansikte.se/ARAGORN/Downloads/aragorn1.2.41.c, and that host answered our request with HTTP 455 on 2026-10-20 (https://github.com/bioconda/bioconda-recipes/blob/master/recipes/aragorn/meta.yaml)
- Copy of the same source: pyaragorn vendors it as vendor/aragorn/aragorn1.2.41.c (pyaragorn README: wraps ARAGORN v1.2.41, "the latest release") (https://github.com/althonos/pyaragorn)
- License GNU GPL v3 or later; maintainer Dean Laslett, e-mail given in the header (aragorn1.2.41.c:6-21)
- Papers: Laslett D, Canback B (2004) Nucleic Acids Res 32:11-16, doi:10.1093/nar/gkh152 (verified by Crossref); the metazoan mitochondrial search is described in the ARWEN paper, doi:10.1093/bioinformatics/btm573 (doi:10.1093/nar/gkh152) (doi:10.1093/bioinformatics/btm573)
- Source header says "ARAGORN (together with ARWEN at last)": the ARWEN code was merged into ARAGORN, which is why `-mt` reproduces ARWEN (see Relevance) (aragorn1.2.41.c:7)

## What it does
- Detects tRNA genes, mitochondrial tRNA genes (`-mt`) and tmRNA genes in FASTA or GenBank nucleotide sequences with heuristic structure/energy scoring (no covariance models).
- With `-mt`/`-mtmam` it uses the ARWEN-derived metazoan mitochondrial search (TV-replacement and D-replacement loops allowed); without them it looks for canonical cloverleaf tRNAs of bacteria/archaea/nuclear genomes and can find C-loop introns (`-i`).
- Reports species from the anticodon through the chosen genetic code; `-w` gives one line per gene.

## Algorithm & key source locations
- Paper description of the standard search: a heuristic that starts from a partially mismatched "GTTC" T-loop hit and builds a cloverleaf around it (T-loop 5-9 nt, T-stem 4-5 bp, TRGYNAA motif 28-85 nt upstream of the T-stem, D-loop 5-11 nt, D-stem 3-4 bp, A-stem 7-9 bp, V-loop 3-25 nt); C-loop introns are not detected unless asked, up to 3000 nt (doi:10.1093/nar/gkh152) (https://academic.oup.com/nar/article/32/1/11/1194039)
- Thresholds in the source: `tRNAthresh` 132.0 (standard), `mtRNAdtthresh` 91.5 (D and T arm present), `mtRNAtthresh` 83.5 (no D-stem), `mtRNAdthresh` 85.0 (no T-stem), `tmRNAthresh` 325.0 (aragorn1.2.41.c:759-763)
- `-e` prints the score as a percentage of the applicable threshold: `100*energy/eref`, with `eref` chosen by whether the D-stem or T-stem is missing (aragorn1.2.41.c:5214-5226)
- `-ps<num>` multiplies the thresholds (default levels 100%, `-ps` alone = 95%); `-rp` flags "possible pseudogenes" (score < 100 or anticodon loop other than 7 nt) and only has an effect when the thresholds are lowered below 100% (aragorn1.2.41.c:12143-12152) (aragorn1.2.41.c:2202-2210)
- Option parsing for metazoan mitochondria: `-mt` sets `mtrna`, selects the composite metazoan genetic code unless a `-gc` was given (`gcfix`); a letter `m` (`-mtmam`) turns the TV-loop search off and selects the vertebrate code; `x` (`-mtx`) clears `mtxdetect`; `d` (`-mtd`) sets `mtcompov` (aragorn1.2.41.c:12324-12357)
- `mtxdetect`: candidates scoring below the threshold are still accepted through "exception" rules (TV-loop, D-loop and cloverleaf branches) unless `-mtx` is given; the help text says `-mtx` = "Same as -mt but low scoring tRNA genes are not reported" (aragorn1.2.41.c:8345) (aragorn1.2.41.c:8797) (aragorn1.2.41.c:9424) (aragorn1.2.41.c:2117-2118)
- `-mt` ignores `-i` and `-sr`, so mitochondrial tRNA genes with introns are not detected (aragorn1.2.41.c:2110-2113)
- Overlap rule (`nearest_trna_gene`): a new candidate that overlaps a stored tRNA by at least `proximity` = 7*MINCTRNALEN/10 = 43 nt (MINCTRNALEN 62) competes with it and only the higher-energy gene is kept; in `-mt` mode this also holds for opposite strands unless `-mtd` is given (aragorn1.2.41.c:6584) (aragorn1.2.41.c:6608-6609) (aragorn1.2.41.c:793)
- Circular sequences: coordinates are wrapped by `sq()`, so a gene across the origin is printed with start > end (aragorn1.2.41.c:6826-6841)
- `-w` species naming: loops of 6 or 8 bases print `mtRNA-???` (or `tRNA-???`) with a 2- or 4-letter anticodon; `-wa` prints `?(Aaa|Bbb)` candidates instead (aragorn1.2.41.c:5826-5860) (aragorn1.2.41.c:2237-2240)
- `-mt` alone also switches the tmRNA search on, because neither `-t` nor `-m` is given; adding `-t` removes it (aragorn1.2.41.c:12418-12421)
- Acceptor end: `-j`, `-jr`, `-jr4` and `-a7` change how the 3' end of the A-stem (NCCA) is predicted or displayed, so reported end coordinates depend on them (aragorn1.2.41.c:2195-2200) (aragorn1.2.41.c:2223)
- Genetic code: `-gc<n>` for GenBank transl_table n, named shortcuts from `-gcstd` to `-gcceph` (Cephalodiscidae), and per-codon changes such as `-gcvert,aga=Trp`; with `-mt` the default is the composite metazoan code, and `-gc` before or after `-mt` gives the same output (tested) (aragorn1.2.41.c:2121-2152) (aragorn1.2.41.c:12325)

## Defaults that matter
- Circular topology (`-c`) and both strands (`-d`) by default; `-l` linear, `-s`/`-sc` one strand (aragorn1.2.41.c:2155-2165)
- Everything detected by default (tRNA and tmRNA); `-t`, `-m` restrict it; `-mt` replaces the tRNA mode with the mitochondrial one (aragorn1.2.41.c:12418-12421) (aragorn1.2.41.c:2102-2109)
- `-ss` (strict spacers) is "Ignored if -mt set"; the non-strict default (3 bp spacer1, 0-2 bp spacer2) "may degrade selectivity" (aragorn1.2.41.c:2191-2194)
- `-tv`: "Do not search for mitochondrial TV replacement loop tRNA genes. Only relevant if -mt used" (aragorn1.2.41.c:2188-2189)
- `-i` default maximum intron length 3000, minimum 0, only for non-mitochondrial searches (aragorn1.2.41.c:2166-2171) (aragorn1.2.41.c:792)
- `-amt<num>`: GenBank input only; length-mismatch report threshold for annotated tRNAs, default 10 bases (aragorn1.2.41.c:2226-2229)
- Input: FASTA or GenBank, results to stdout; the internal sequence-window constant LSEQ is 20000 (aragorn1.2.41.c:850) (aragorn1.2.41.c:2099-2101)

## Known bugs / issues
- No public issue tracker for ARAGORN itself was found (search not exhaustive 未确证); the header asks for bug reports by e-mail and thanks several users for finding bugs (aragorn1.2.41.c:23-34)
- The 2004 paper evaluated archaeal, bacterial and nuclear tRNAs and chloroplast tmRNAs; it does not mention mitochondrial genomes, so ARAGORN `-mt` performance claims come from the ARWEN paper and the MiTFi comparison (doi:10.1093/nar/gkh152) (doi:10.1093/nar/gkr1131)
- Random-sequence false-positive rate of the standard mode: 1.7 per Gb at 40% G+C, 14.3 per Gb at 60%, 396.9 per Gb at 80% (not measured for `-mt`) (https://academic.oup.com/nar/article/32/1/11/1194039)
- `-mt` is the same engine as ARWEN: on six test mitogenomes (one annelid, four anthozoans, one sponge; a seventh record, a 1.5 kb COX1 fragment, gave no calls in any tool) the species and coordinates printed by `aragorn -mt -w` were identical to `arwen -w`, and `aragorn -mtx -w` identical to `arwen -mtx -w` (observed 2026-10-20) (aragorn1.2.41.c:7) (aragorn1.2.41.c:481-483)
- Strictness matters: on NC_000931 `-mt` gave 27 calls, `-mtx` 25, `-mt -mtd` 31 (four extra opposite-strand overlaps); on NC_000933 (Metridium) `-mt` 6 vs `-mtx` 5; on NC_010217 10 vs 9 (observed 2026-10-20) (aragorn1.2.41.c:2117-2120) (aragorn1.2.41.c:6608-6609)
- The `-e` score does not separate true from false calls: NC_000931 true calls 100.2-115.3, false calls 94.6-102.1; NC_006902 true 103.6 and 115.8, false 100.7-105.5 (observed 2026-10-20) (aragorn1.2.41.c:5214-5226)
- Calls with a 6- or 8-base anticodon loop have no species in `-w` output (`???`), so a downstream parser that needs an amino acid silently loses them; e.g. NC_006902 position 17607-17664 (observed 2026-10-20) (aragorn1.2.41.c:5845-5855)
- ARAGORN reports tmRNA genes too if only `-mt` is given (the run prints "Searching for tmRNA genes"); in animal mtDNA a tmRNA call is expected to be false (observed 2026-10-20) (aragorn1.2.41.c:12418-12421)

## Limitations by taxon
- Metazoan mitochondria with many lost tRNAs (Cnidaria, sponges, placozoans): in our benchmark ARAGORN gave 462 false positives against 183 true tRNAs on 85 RefSeq mitogenomes (91.8% sensitivity); annelids 117 false positives on 40 mitogenomes (96.5%) (README.md:670)
- Test genomes with 2 RefSeq tRNAs: `aragorn -mt` reported 6 (Metridium senile), 5 (Montipora cactus), 5 (Nematostella sp.), 8 (Exaiptasia diaphana) and 10 for Plakinastrella (6 RefSeq tRNAs); 22 RefSeq tRNAs in Platynereis dumerilii gave 27 calls, 21 of them matching RefSeq (observed 2026-10-20; overcalling follows from the exception rules) (aragorn1.2.41.c:8345)
- Arachnida and tRNAs lacking the 3' end: the ARWEN engine recovers more of them than tRNAscan-SE but with many false calls (see arwen.md; evaluation of ARWEN: sensitivity 0.977, precision 0.869 on 1,876 RefSeq mitogenomes) (doi:10.1093/nar/gkr1131)
- Nematodes with armless tRNAs: no ARAGORN-specific data found in the sources read 未确证 (doi:10.1093/nar/gkr1131)
- Non-metazoan mitochondria (fungi, plants, protists): `-mt` is documented as "Search for Metazoan mitochondrial tRNA genes"; its suitability for other lineages is not documented 未确证 (aragorn1.2.41.c:2110-2113)

## Relevance to our pipeline
- `mito_trna.py` runs `aragorn -mt -gc<code> -w` (mito_trna.py:127-128) and ARWEN as `arwen -gc<code> -mtx -w` (mito_trna.py:123-124). Both are the same engine, so their agreement is not independent evidence (the pipeline already counts them as one family, mito_trna.py:248-252). The only difference is `-mt` (exception rules on) versus `-mtx` (off): this very probably explains why ARAGORN gives more false calls than ARWEN in the README table (462 vs 348 on cnidarians, 117 vs 79 on annelids); on our six test mitogenomes the outputs were identical whenever the flags matched (not re-checked on the 125 benchmark genomes). For a consistent second-method set, pass `-mtx` to ARAGORN too, or drop one of the two.
- Genetic code: pass the same `--genetic_code` as for MITOS2 (`-gc5`, `-gc4`); with `-mt` and `-gc` the order of the two flags does not matter.
- Parser: `LASLETT` requires a 3-letter anticodon and an amino-acid name; `mtRNA-???` calls (6/8 anticodon loops) are dropped (mito_trna.py:96-109). That is acceptable for supporting evidence but means ARAGORN cannot support tRNAs with atypical anticodon loops.
- Origin-spanning genes are written with start > end; the parser converts them to extended coordinates (mito_trna.py:101-118), consistent with `sq()` wrapping in the source.
- `-mt -t` avoids the useless tmRNA search (runtime difference was negligible, 0.83 s vs 0.85 s for 15.6 kb).
- Boundaries: in our NC_000931 test only 10 of 21 matching calls had boundaries identical to RefSeq (MITOS2 22 of 22), so keep taking boundaries from MITOS2 when it found the locus.
- Score filter: `-e` scores overlap between true and false calls in the test genomes, so a score cut-off on ARAGORN output is not a safe replacement for the MITOS2 requirement.

## References
- Laslett D, Canback B (2004) ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences. Nucleic Acids Res 32:11-16 (doi:10.1093/nar/gkh152)
- Laslett D, Canback B (2008) ARWEN: a program to detect tRNA genes in metazoan mitochondrial nucleotide sequences. Bioinformatics 24(2):172-175 (doi:10.1093/bioinformatics/btm573)
- Jühling F et al. (2012) Improved systematic tRNA gene annotation ... (MiTFi). Nucleic Acids Res 40(7):2833-2845 (doi:10.1093/nar/gkr1131)
- ARAGORN 1.2.41 source (https://github.com/althonos/pyaragorn/blob/main/vendor/aragorn/aragorn1.2.41.c)
- Bioconda recipe (https://github.com/bioconda/bioconda-recipes/blob/master/recipes/aragorn/meta.yaml)
