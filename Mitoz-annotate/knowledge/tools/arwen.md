# ARWEN

Every bullet outside "What it does" / "Relevance to our pipeline" ends with a source:
`(path/in/source.py:123)`, `(doi:10.xxxx/...)`, `(https://...)`, or is marked 未确证.
Check: `python3 benchmark/check_knowledge.py knowledge/`

Path convention: `arwen1.2.3.c:N` = line N of the source file `arwen1.2.3.c` (github.com/RemiAllio/MitoFinder, directory `arwen/`, md5 2f43d6cbd77123f4e224abb04a56d815 as downloaded 2026-10-20); `aragorn1.2.41.c:N` = line N of the ARAGORN 1.2.41 source (see aragorn.md).

## Version checked
- ARWEN 1.2.3 (source header "Version 1.2.3, May 19th, 2008"; the program banner and `-h` print only "ARWEN v1.2"), checked 2026-10-20 (arwen1.2.3.c:4) (arwen1.2.3.c:6645)
- Local binary /root/miniconda3/envs/mitotrna/bin/arwen (175,016 bytes, file date Oct 10 12:10) prints "ARWEN v1.2" and contains the string `arwen1.2.3.c`; there is no conda package record for it, so who compiled it and with which flags is 未确证 (`arwen -h` banner)
- ARWEN is not on bioconda; the source is distributed inside MitoFinder and MitoPilot, which compile it with `gcc` (MitoPilot: `-O2 -ffast-math -lm`) (https://github.com/RemiAllio/MitoFinder/blob/master/arwen/arwen1.2.3.c) (https://github.com/NMNH-OceanDNA/MitoPilot/blob/main/inst/native/install_mitopilot_native.sh)
- License: GNU GPL version 2; authors Dean Laslett and Bjorn Canback (arwen1.2.3.c:4-30)
- Paper: Laslett D, Canback B (2008) Bioinformatics 24(2):172-175, doi:10.1093/bioinformatics/btm573 (verified by Crossref: published online 2007-11-22) (doi:10.1093/bioinformatics/btm573)

## What it does
- Heuristic (not covariance-model) detector of metazoan mitochondrial tRNA genes in FASTA nucleotide sequences, including mitochondrial tRNAs with a TV-replacement loop (no T-arm) or a D-replacement loop (no D-arm).
- Prints genes in text form: tRNA species (from the anticodon and the chosen genetic code), start/end, strand ("c" = complementary), anticodon position and anticodon triplet.
- ARAGORN 1.2.x contains the same mitochondrial search routines behind `-mt`; ARWEN is the earlier stand-alone program.

## Algorithm & key source locations
- Search strategy: a heuristic hunt for a C-arm (anticodon) hairpin with a 5-6 bp stem and 6-8 base loop; upstream a D-arm (2-5 bp stem, 3-17 base loop), downstream a T-arm (2-7 bp stem, 3-31 base loop), and an A-stem (5-8 bp) at both ends; a candidate is accepted when at least three of the four structures are present (https://academic.oup.com/bioinformatics/article/24/2/172/228888)
- Three scoring algorithms: standard cloverleaf, D-replacement loop and TV-replacement loop; cutoffs lowest for D-replacement, then TV-replacement, then cloverleaf; bonus for a TA spacer, a "CT<nnn>AA" C-arm and GTTC homology in the T-arm; GC pairs not penalised, AT/GT/non-bonding pairs penalised in increasing order (https://academic.oup.com/bioinformatics/article/24/2/172/228888)
- Score thresholds in the source: `mtRNAdtthresh` 91.5 (cloverleaf with D and T arm), `mtRNAtthresh` 83.5 and `mtRNAdthresh` 85.0 (one-arm variants) (arwen1.2.3.c:481-483)
- `-mtx` switches off the "exception" rules that rescue candidates scoring below the thresholds (the three `if (!mtxdetect)` branches for TV-loop, D-loop and full cloverleaf candidates) (arwen1.2.3.c:4408) (arwen1.2.3.c:4859) (arwen1.2.3.c:5485) (arwen1.2.3.c:6923)
- Length constants: MINCTRNALEN 62, MAXCTRNALEN 110, A-stem extension 9, MAXINTRONLEN 3000 (introns are not searched in ARWEN, unlike ARAGORN `-i`) (arwen1.2.3.c:506) (arwen1.2.3.c:513-518) (https://academic.oup.com/bioinformatics/article/24/2/172/228888)
- Overlap rule: a new candidate is compared with stored genes; if it overlaps a stored tRNA by at least `proximity` = 7*62/10 = 43 nt (integer arithmetic) the candidate with the lower energy is dropped, and in mitochondrial mode this applies to overlaps between opposite strands too (ARWEN has no `-mtd` switch to allow them) (arwen1.2.3.c:2781-2790) (arwen1.2.3.c:2925-2945)
- Species naming comes from the anticodon through the chosen genetic code; for 6- or 8-base anticodon loops (`cloop` 6 or 8) `-w` prints the species as `mtRNA-???` and an anticodon of only 2 or 4 bases (arwen1.2.3.c:2668-2690)
- Genetic code: `-gc<n>` selects a GenBank transl_table; the composite metazoan code (`-gcmet`) reports more than one possible iso-acceptor where metazoan codes differ (arwen1.2.3.c:6668-6686) (https://academic.oup.com/bioinformatics/article/24/2/172/228888)
- Topology and strands: circular (`-c`) and both strands (`-d`) are the defaults; `-tv` disables the TV-replacement search; `-rp` reports low-scoring genes as possible pseudogenes (arwen1.2.3.c:6694-6708)
- Acceptor end: `-j` / `-jr` / `-jr4` control how the 3' end of the A-stem (NCCA) is displayed or allowed to diverge, so reported end coordinates depend on these switches and on the predicted acceptor length (arwen1.2.3.c:6715-6720)

## Defaults that matter
- Search both strands, circular topology, TV-replacement loops searched, standard threshold level (no `-ps`-like switch exists in ARWEN) (arwen1.2.3.c:6694-6708)
- `-mtx` = low scoring tRNA genes are not reported (arwen1.2.3.c:6666)
- Input must be FASTA; several sequences are allowed; results go to STDOUT (arwen1.2.3.c:6650-6657)
- Without `-w` the program prints a header ("Assuming circular topology, search wraps around ends", "Searching both strands", the genetic code in use) and one numbered entry per gene; `-w` gives one line per gene with species, `[start,end]`, anticodon position and `(anticodon)` (observed 2026-10-20 with the local binary) (arwen1.2.3.c:2668-2700)
- Internal buffers are sized 2*LSEQ+WRAP with LSEQ = 20000 (arwen1.2.3.c:570) (arwen1.2.3.c:6340)

## Known bugs / issues
- Designed for sensitivity: "ARWEN buys its increased sensitivity at the expense of a substantial false discovery rate"; the authors themselves state it is "not suitable for use on longer cytoplasmic genomes" (doi:10.1093/nar/gkr1131) (https://academic.oup.com/bioinformatics/article/24/2/172/228888)
- Original paper numbers: on 23 RefSeq mitogenomes 99.4% sensitivity but 115 false positives (tRNAscan-SE 72.1% and 8); raising the thresholds to 8 false positives dropped sensitivity to 81.7%; "sixty-four percent of the falsely predicted tRNAs are derived from 6 genomes"; random-sequence calibration about 186.6 false positives per Mb, "almost 3 false positives per genome on average" (https://academic.oup.com/bioinformatics/article/24/2/172/228888)
- Independent evaluation on 1,876 RefSeq metazoan mitogenomes (40,521 annotated tRNAs): ARWEN sensitivity 0.977, precision 0.869 (39,569 TP, 5,957 FP, 952 FN); MiTFi 0.986 and 0.979; the genomes came from NCBI RefSeq version 39 (doi:10.1093/nar/gkr1131)
- Output ambiguity: tRNAs with 6- or 8-base anticodon loops are listed as `mtRNA-???` with a 2- or 4-letter anticodon; the species cannot be read from this output (use ARAGORN `-wa` for candidate species) (arwen1.2.3.c:2679-2683) (aragorn1.2.41.c:2237-2240)
- Boundaries are not RefSeq-like: the gene length `nbase` includes the A-stem plus `aext` extra acceptor-tail bases predicted from the 3' sequence, and the stop is start + nbase - 1; in our run on NC_000931 (21 calls matching one of 22 RefSeq tRNAs) the 5' end was identical in 11 and 1 nt upstream in 10, and the 3' end identical in 10, 1 nt downstream in 10 and 2 nt downstream in 1, never shorter (observed 2026-10-20) (arwen1.2.3.c:6044-6046) (arwen1.2.3.c:2934)
- No separate issue tracker is known for ARWEN; bugs are reported to the authors by e-mail (arwen1.2.3.c:10-14)

## Limitations by taxon
- Taxa with few tRNA genes (Cnidaria, Chaetognatha, sponges, ctenophores): the false-positive burden is highest because the program always tries to explain every stem-loop; in our benchmark ARWEN gave 348 false positives against 183 true tRNAs on 85 cnidarian/sponge/placozoan RefSeq mitogenomes (README.md:670)
- Annelida: 79 false positives on 40 RefSeq mitogenomes (880 tRNAs), 95.7% sensitivity (README.md:670)
- Arachnida: ARWEN found 82.8% of annotated tRNAs (MiTFi 89.7%, tRNAscan-SE 50.4%); for the two families examined it made 233 false-positive calls (MiTFi 62, tRNAscan-SE 22) (doi:10.1093/nar/gkr1131)
- Truncated tRNAs lacking the 3' end (Cecidomyiidae): the MiTFi paper quotes "ARWEN reported 21 and tRNAscan-SE recovered only 7" (denominator not read) (doi:10.1093/nar/gkr1131)
- Nematodes (Enoplea: many tRNAs lack the T-domain, some lack D and T): ARWEN's TV-replacement and D-replacement scoring covers a lost single arm, but no ARWEN-specific measurement for armless nematode tRNAs was found in the sources read 未确证 (doi:10.1093/nar/gkr1131)
- Mammalian mitochondria: `-mtmam` fixes the TV-loop search off and uses the vertebrate code; the paper reports 100% sensitivity on 677 Mamit-tRNA sequences (arwen1.2.3.c:6915-6921) (https://academic.oup.com/bioinformatics/article/24/2/172/228888)

## Relevance to our pipeline
- `mito_trna.py` calls `arwen -gc<code> -mtx -w <fasta>` and parses `mtRNA-Xxx c[start,end] pos (anticodon)` with `LASLETT` (mito_trna.py:123-124, 96-109). `-mtx` is the stricter setting; ARAGORN is called with `-mt` (no x), so the two "same family" tools are run at different strictness.
- Because ARWEN is an overlap-tolerant heuristic with high false-positive rate in tRNA-poor taxa, the pipeline only uses it as supporting evidence: a MITOS2 call is kept alone, an ARWEN-only call needs a second independent family (ARWEN and ARAGORN count as one family) (README.md:640-645, 676).
- Our test run 2026-10-20 (`arwen -gc5 -mtx -w`, NC_000931 Platynereis dumerilii, 22 RefSeq tRNAs): 25 calls, 21 overlapping a RefSeq tRNA on the same strand, only 10 with identical start and end (the others 1 nt longer at the 5' end and/or 1-2 nt longer at the 3' end), 4 false positives (Ser tct, Pro agg, Leu aag, Tyr ata), versus MITOS2 22 of 22 with identical boundaries. On NC_006902 (Montipora cactus, code 4, 2 RefSeq tRNAs): 5 calls, trnM and trnW with identical boundaries, plus Ser(tga), Ile(aat) and one `???` call (2-base anticodon "tg").
- Do not take ARWEN boundaries over MITOS2 boundaries; the pipeline already takes boundaries from MITOS2 first.
- `???` rows: the regex requires a 3-letter anticodon, so 6/8-loop calls are dropped, which is acceptable for supporting evidence but means ARWEN contributes nothing for tRNAs with atypical anticodon loops.

## References
- Laslett D, Canback B (2008) ARWEN: a program to detect tRNA genes in metazoan mitochondrial nucleotide sequences. Bioinformatics 24(2):172-175 (doi:10.1093/bioinformatics/btm573)
- Jühling F et al. (2012) Improved systematic tRNA gene annotation allows new insights into the evolution of mitochondrial tRNA structures and into the mechanisms of mitochondrial genome rearrangements. Nucleic Acids Res 40(7):2833-2845 (doi:10.1093/nar/gkr1131)
- Laslett D, Canback B (2004) ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences. Nucleic Acids Res 32:11-16 (doi:10.1093/nar/gkh152)
- ARWEN source arwen1.2.3.c (https://github.com/RemiAllio/MitoFinder/blob/master/arwen/arwen1.2.3.c)
