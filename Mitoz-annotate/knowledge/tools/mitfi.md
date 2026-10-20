# MiTFi (mitochondrial tRNA finder, as bundled in MitoZ)

Every bullet outside "What it does" / "Relevance to our pipeline" ends with a source:
`(path/in/source.py:123)`, `(doi:10.xxxx/...)`, `(https://...)`, or is marked 未确证.
Check: `python3 benchmark/check_knowledge.py knowledge/`

## Version checked
- MiTFi v0.1 (banner printed by the jar), the copy bundled in MitoZ 3.6 under mitoz/annotate/script/mitfi/, checked 2026-10-20 (mitoz/annotate/script/mitfi/README.TXT:1)
- Installed jar md5 62c373a654aded66c3e5c5459f67eef5 equals the md5 the MitoZ author printed for the MitoZ 2.3 jar in 2019, so the bundled jar is unchanged since then; jar entries are dated 2011-04-27 (https://github.com/linzhi2013/MitoZ/issues/20)
- Search engine inside the bundle is Infernal 1.0.2 (`infernal-1.0.2.binaries/cmsearch.{linux,darwin_i386,darwin_arm64}`, linked as `cmsearch` by detect_platform.sh); MitoZ rRNA search uses a different, newer Infernal 1.1.1 (mitoz/annotate/script/mitfi/detect_platform.sh:19)
- The installed bundle holds only the compiled jar (2 classes: `mitfi.Main`, `mitfi.tRNA`, no .java), the 22 CMs, gc.prt, config and Infernal binaries; upstream paper is Jühling et al. 2012 (doi:10.1093/nar/gkr1131). The upstream MiTFi download page could not be fetched (certificate error on the legacy host), so upstream version history is 未确证

## What it does
- Finds the 22 metazoan mitochondrial tRNA genes in a mitogenome FASTA by running Infernal covariance-model searches (one CM per amino acid, two each for Leu and Ser), predicting an anticodon for every hit, and choosing among overlapping hits by E-value. MitoZ runs it once per input sequence and turns the table into feature-table lines.

## Algorithm & key source locations
- Pipeline described in the paper: all 22 CMs are run against the genome; for every locus only the hit of the CM with the smallest E-value is kept; "no score cutoffs are used for the 22 top hits at this point"; predictions overlapping by at most 10 nt count as distinct loci (doi:10.1093/nar/gkr1131)
- Second step in the paper: copies of tRNA genes are searched in the remaining locations with an E-value cutoff of 0.001, chosen because true positives are nearly unaffected at E = 0.001 while false positives drop to a constant level (doi:10.1093/nar/gkr1131)
- Anticodon prediction: count interior stem-loops; if only two, scan loops for 7 unpaired nt (anticodon loop), else take the loop closest to the sequence mean, and consider loop size 9 if no 7-nt loop exists; if no anticodon can be assigned the hit is labelled with the amino acid of the CM that found it (doi:10.1093/nar/gkr1131)
- Amino acid label comes from anticodon plus the chosen NCBI genetic code (`-code`, table read from gc.prt), so the label can differ from the CM that matched; example from a real run: anticodon UAG labelled `L1` from model `Metazoa_L2.cm` (https://github.com/linzhi2013/MitoZ/issues/208)
- CM files: 22 files `Metazoa_{A,C,D,E,F,G,H,I,K,L1,L2,M,N,P,Q,R,S1,S2,T,V,W,Y}.cm`, built with Infernal 1.0 format; seed sizes 33-69 sequences in the paper, e.g. Metazoa_A.cm NSEQ 36 CLEN 65, Metazoa_L2.cm NSEQ 45 CLEN 66 (mitoz/annotate/script/mitfi/models/Metazoa_A.cm:1)
- Jar entry points seen with `javap`: `Main.findAnticodon`, `findAnticodon_only7`, `InfernalParser`, `addHits`, `reduceHits`, `overlapping`; `tRNA` objects carry fields `intron`, `istart`, `istop`, `astart`, `astop`, `score_relative`; what the intron fields do is 未确证 (mitoz/annotate/script/mitfi/mitfi.jar)
- Command line as printed by `java -jar mitfi.jar -h` (observed 2026-10-20): `-fasta -structure -hits -progress -bstrands -overlap n[10] -top -cm f -code c[5] -codefile f -onlycutoff -evalue n[0.001] -cores m`; `-onlycutoff` is described there as "hits are selected by cutoff only, no top-22 round" (mitoz/annotate/script/mitfi/README.TXT:24-25)
- Tab output columns: `#header start stop score evalue AC AA model strand`; AC is the anticodon in RNA letters (U not T); start < stop even on the minus strand, strand given in last column (https://github.com/linzhi2013/MitoZ/issues/208)
- How MitoZ converts it: one temp FASTA per record, `cd mitfi; ./detect_platform.sh; java -Xmx2048m -jar mitfi.jar -cores 1 -code <table> -evalue 0.001 -onlycutoff <file>` appended to `<fasta>.trna` (mitoz/annotate/annotation.py:257-279)
- Conversion script: removes `1`/`2` from the amino-acid label (so L1/L2 become L, S1/S2 become S), lower-cases the anticodon, writes gene name `trn<AA>(<anticodon>)` and product `tRNA-<Three-letter>`; minus-strand hits are written with start/stop swapped (mitoz/annotate/script/tRNA_ft_mitfi.py:58-73)
- Single core is forced because multi-core MiTFi returned no results (author's code comment) (mitoz/annotate/annotation.py:264)

## Defaults that matter
- Config file: `evalue = 0.001`, `overlap = 10`, `code = 5`, `cores = 1`, `infernalCall = ./cmsearch` (mitoz/annotate/script/mitfi/mitfi_config.txt:8-24)
- MitoZ passes `-code` = the annotate `--genetic_code` (auto: Chordata 2, Echinodermata 9, Porifera 4, all other listed clades 5), `-evalue 0.001`, `-onlycutoff`, `-cores 1`, heap 2 GB (mitoz/annotate/annotation.py:274-278)
- MitoZ paper states it outputs tRNAs with E <= 0.001 using `-cores 1 -evalue 0.001 -onlycutoff -code 2/5` (doi:10.1093/nar/gkz173)
- Inference, not tested: with `-onlycutoff` the first "best hit per CM, no cutoff" round of the paper is skipped, so a tRNA whose best hit has E > 0.001 is never reported; this matches reports of missing tRNAs but the causal link is 未确证 (https://github.com/linzhi2013/MitoZ/issues/74)
- Run time: MiTFi takes about 80 % of total `mitoz annotate` time on a 15 kb annelid genome (137 of 168 s) (temp_scripts/Mitoz-annotate/scripts/mitoz_fast.py:4-7)

## Known bugs / issues
- Java `StringIndexOutOfBoundsException ... mitfi.Main.main(Main.java:221)` on a CentOS 7 host with Zulu OpenJDK 11, unresolved; the reporter switched to Singularity (https://github.com/linzhi2013/MitoZ/issues/20)
- Bundled `cmsearch.linux` fails to execute on some hosts ("Exec format error", "Bad address", WSL), so MiTFi finds no tRNAs, and MitoZ continues without error; the same run also doubled/split the s-rRNA in issue 187; fix is to replace the binary or compile Infernal 1.0.2 (https://github.com/linzhi2013/MitoZ/issues/187)
- Same "no tRNA in the result" symptom on Windows Subsystem for Linux and CentOS 8 conda installs, solved by rebuilding `cmsearch` (https://github.com/linzhi2013/MitoZ/issues/208)
- Fix for the broken-binary cases is shipped with the bundle: compile Infernal 1.0.2 and use its `cmsearch` in place of `cmsearch.linux` (mitoz/annotate/script/mitfi/mgl_readme.txt:1-12)
- Missing tRNAs for spiders and nematodes (about 10 absent in >14 kb genomes) while MITOS finds them; maintainer: the CM set has not been updated for years, build your own CMs (https://github.com/linzhi2013/MitoZ/issues/74)
- Marsupial tRNA-Asp with GCC anticodon is reported as a second tRNA-Gly; no option to edit expected anticodons (https://github.com/linzhi2013/MitoZ/issues/237)
- User reported tRNA missing in MitoZ but found by tRNAscan-SE, solved by replacing `cmsearch.linux` (https://github.com/linzhi2013/MitoZ/issues/225)
- Leu/Ser copies share names after the digit is stripped; they differ only by anticodon in the gene name, and downstream code that keys on `trnL`/`trnS` alone will merge them (mitoz/annotate/script/tRNA_ft_mitfi.py:58)
- tRNA vs CDS overlap is never reconciled in MitoZ: an ND1 call covered a whole trnS in one report; maintainer pointed to MiTFi for tRNA calls (https://github.com/linzhi2013/MitoZ/issues/146)

## Limitations by taxon
- Paper: frequent tRNA gene loss in basal Metazoa and loss of individual arms (D or T arm) in Arthropoda; CMs include "some of the aberrant tRNA structures" but cannot cover all (doi:10.1093/nar/gkr1131)
- NCBI expects 22 tRNAs but lists Cnidaria and Porifera as known exceptions, so a "missing" tRNA there is not necessarily an error (https://www.ncbi.nlm.nih.gov/genbank/organelle_submit/)
- Reported weak taxa: spiders, nematodes (issue 74), marsupials (issue 237); annelid-specific performance is 未确证
- MitoZ control-region detection demands 22 distinct tRNA names, so any family with a lost tRNA gets no control-region feature (mitoz/annotate/script/find_contro_region.py:118-176)

## Relevance to our pipeline
- MiTFi is the slowest step of `mitoz annotate`; `scripts/mitoz_fast.py` swaps it for MITOS2 (same Metazoa CMs) and writes MiTFi's own `.trna` format so `tRNA_ft_mitfi.py` and the rest run unchanged; keep the column order `start stop score evalue AC AA model strand` and RNA-letter anticodons when writing replacement results.
- Because MitoZ cuts at E <= 0.001 without the top-22 rescue, absent tRNAs in worm genomes are expected; `mito_trna.py` (multi-tool, MITOS2 first) is the right place to recover them, and a tRNA count below 22 should not be forced.
- Gene names are `trnX(anticodon)` with lower-case RNA letters (e.g. `trnL(uaa)`), unlike NCBI-style DNA anticodons; any join with GenBank references must normalise this.
- tRNA boundaries from MiTFi are kept by `mito_curate.py` (criterion v), so MiTFi/MITOS2 boundary choices directly determine where T/TA stop codons can be placed.

## References
- Jühling F, Pütz J, Bernt M, Donath A, Middendorf M, Florentz C, Stadler PF (2012) Improved systematic tRNA gene annotation allows new insights into the evolution of mitochondrial tRNA structures and into the mechanisms of mitochondrial genome rearrangements. Nucleic Acids Res 40(7):2833-2845 (doi:10.1093/nar/gkr1131)
- Meng G, Li Y, Yang C, Liu S (2019) MitoZ: a toolkit for animal mitochondrial genome assembly, annotation and visualization. Nucleic Acids Res 47(11):e63 (doi:10.1093/nar/gkz173)
- MitoZ issue tracker (https://github.com/linzhi2013/MitoZ/issues)
