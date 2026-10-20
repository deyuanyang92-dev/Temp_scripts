# MitoFinder (Allio et al.)

Every bullet outside "What it does" / "Relevance to our pipeline" ends with a source:
`(path/in/source.py:123)`, `(doi:10.xxxx/...)`, `(https://...)`, or is marked 未确证.
Check: `python3 benchmark/check_knowledge.py knowledge/`

Path prefix `MitoFinder/` = https://github.com/RemiAllio/MitoFinder at commit e8ab3ad (shallow clone, read only; nothing was installed or run).
`Mitoz-annotate/` = our pipeline directory.

## Version checked
- MitoFinder v1.4.2 (release 2024-04-11; repository HEAD e8ab3ad of 2025-09-02, message "Update README.md"; the main script and `geneChecker_fasta.py` say version 1.4.2); checked 2026-10-20 (https://github.com/RemiAllio/MitoFinder/releases) (MitoFinder/mitofinder:1-2)
- Paper: Allio R, Schomaker-Bastos A, Romiguier J, Prosdocimi F, Nabholz B, Delsuc F. 2020, Mol Ecol Resour 20:892-905 (doi:10.1111/1755-0998.13160)
- Not run in this session: no Python 2 on this machine and `conda search -c bioconda mitofinder` did not finish within 110 s, so package facts come from the anaconda.org API and the bioconda recipe (https://api.anaconda.org/package/bioconda/mitofinder)

## What it does
- Pipeline that assembles a mitogenome from trimmed reads (MEGAHIT, metaSPAdes or IDBA-UD), finds the mitochondrial contig(s) by BLAST against user-supplied GenBank mitogenome reference(s), and annotates protein-coding genes, rRNAs and tRNAs.
- Can also skip assembly and annotate contigs of an existing assembly (`-a`), which is the mode of interest for benchmarking against a finished mitogenome FASTA.
- Annotation is reference-based: gene boundaries come from BLAST hits to the reference genes, followed by a short start/stop-codon search; tRNAs come from MiTFi (default), ARWEN or tRNAscan-SE.
- Outputs per contig: GenBank, GFF3, tbl (for NCBI submission), gene FASTA (nucleotide and amino acid), PNG map.

## Algorithm & key source locations
- Reference input is a GenBank file; for each CDS and rRNA the `gene` (else `product`) qualifier is the name; names are normalised through an alias list (coi, cox1, "cytochrome c oxidase subunit I", ...) and must end up as one of COX1, COX2, COX3, CYTB, ND1-6, ND4L, ATP6, ATP8, rrnL, rrnS, otherwise the run stops unless `--ignore` or `--new-genes` is given (MitoFinder/mitofinder:520-560) (MitoFinder/mitofinder:638) (MitoFinder/mitofinder:713-716)
- Contig identification: `blastn` of the reference gene sequences against the assembly with `-evalue 1e-5` and `-perc_identity 50`; a contig is kept only if its length is between `--min-contig-size` (1000) and `--max-contig-size` (25000) (MitoFinder/mitofinder:848) (MitoFinder/mitofinder:869) (MitoFinder/mitofinder:84-85) (MitoFinder/mitofinder:100-103)
- Protein-coding genes: `blastx -query_gencode <table> -seg no` of the contig against the reference protein set; an HSP is accepted when amino-acid identity >= `--blast-identity-prot` (40) and its length is >= `--blast-size` (30 %) of the reference gene (MitoFinder/geneChecker_fasta.py:131) (MitoFinder/geneChecker_fasta.py:147-150)
- The paper (version 1.0) describes tblastx for CDS, blastn for 12S/16S and ARWEN 1.2 for tRNAs, validated on ants against 29 GenBank ant mitogenomes; the current code uses blastx (https://pmc.ncbi.nlm.nih.gov/articles/PMC7497042/) (MitoFinder/geneChecker_fasta.py:131)
- rRNA: `blastn -word_size 8 -perc_identity 50 -max_hsps 5 -dust no` against the reference rRNAs (MitoFinder/geneChecker_fasta.py:216-229)
- Start/stop refinement: after the BLAST hit, the end is moved codon by codon (up to `--nwalk` 5 codons = 30 nt, both directions) until a start codon (5') and a stop codon (3') of the chosen NCBI table is found; strand -1 genes are handled by the mirrored code (MitoFinder/genbankOutput.py:147-300) (MitoFinder/mitofinder:86-87)
- Incomplete stop codons (T / TA completed by polyadenylation) have no special handling in the walk, which only looks for complete stop codons of the table; effect on animal genes 未确证 (MitoFinder/genbankOutput.py:213-240)
- tRNA step: `-t mitfi` (default, `java -jar mitfi.jar -code <table>`, bundled Infernal 1.0.2), `-t arwen`, or `-t trnascan` which runs `tRNAscan-SE -X 7 [-g <code file>] -O` (organelle mode, cutoff fixed to 7 bits) (MitoFinder/mitofinder:69) (MitoFinder/mitofinder:348) (MitoFinder/tRNAscanChecker.py:262-320)
- Introns: only with `--allow-intron` (and `--intron-size`, default 5000): BLASTX pieces of one gene within the allowed distance are merged (`geneChecker_fasta_gaps.py`); the README states that start/stop codons are searched but "there is no search for intronic boundaries" and that a close reference is required (MitoFinder/README.md:479-482) (MitoFinder/geneChecker_fasta_gaps.py:150) (MitoFinder/mitofinder:1083-1084)
- Circularisation check: blast of the two contig ends (`--circular-size` 45, `--circular-offset` 200); informational, run only when one contig is found (MitoFinder/mitofinder:963) (https://github.com/RemiAllio/MitoFinder/issues/62)
- No rotation of the sequence to a standard start: the code that reordered the record to start at tRNA-Phe is inside a commented-out block; `--adjust-direction` only reverse-complements to match the reference (MitoFinder/geneChecker_fasta.py:529-554) (MitoFinder/mitofinder:89)
- Bundled third-party code: BLAST+ 2.2.29, an old Python 2 Biopython copy, Infernal 1.0.2 (MiTFi) and 1.1.3 + tRNAscan-SE 2.0.6 compiled by `install.sh`, ARWEN source (MitoFinder/blast/ncbi_package_info:2) (MitoFinder/install.sh:1-50)

## Defaults that matter
- Genetic code `-o` has no default and must be given (NCBI tables 1-25 accepted) (MitoFinder/mitofinder:210-212)
- `-e` BLAST E-value 1e-5; `--blast-identity-nucl` 50; `--blast-identity-prot` 40; `--blast-size` 30; `--nwalk` 5; `--min-contig-size` 1000; `--max-contig-size` 25000; `--intron-size` 5000; `-t mitfi` (MitoFinder/mitofinder:69-110)
- A genome longer than 25 kb is skipped unless `--max-contig-size` is raised (issue 66 uses a huge value for a long contig) (https://github.com/RemiAllio/MitoFinder/issues/66) (MitoFinder/mitofinder:869)
- For an existing assembly (`-a`) `--rename-contig` defaults to "no"; for de novo it is "yes" (MitoFinder/mitofinder:258-261)
- `--new-genes` is needed for non-standard genes (e.g. fungal rps3); all references must use the same names for them (MitoFinder/mitofinder:91)
- Annotation-only command: `mitofinder -j ID -a assembly.fasta -r reference.gb -o <table> -p <threads> -m <GB>` (MitoFinder/README.md:288)
- Restart behaviour: the same command skips finished steps; `--override` recomputes (MitoFinder/README.md:292-294)
- Requirements stated by the authors: Linux or macOS, automake, autoconf, gcc, Java for MiTFi, Python 2.7; Singularity image recommended (MitoFinder/README.md:17-18) (MitoFinder/README.md:44) (MitoFinder/README.md:93-102)

## Known bugs / issues
- Python 2 only: `mitofinder -h` fails with a syntax error on the `print` statement if `python` is Python 3; fix is `#!/usr/bin/python2.7` in every script (https://github.com/RemiAllio/MitoFinder/issues/54) (MitoFinder/mitofinder:1)
- Python 3 ports exist only as third-party forks (jacquefm, Adamtaranto); the upstream refactor request is still open (https://github.com/RemiAllio/MitoFinder/issues/19) (https://github.com/RemiAllio/MitoFinder/issues/55)
- Bioconda: anaconda.org lists only versions 1.4 and 1.4.1 (linux-64 and osx-64, all `py27`); the recipe on master says 1.4.2 with `skip: True # [not py2k]`; why no 1.4.2 package appears is 未确证 (https://api.anaconda.org/package/bioconda/mitofinder) (https://github.com/bioconda/bioconda-recipes/blob/master/recipes/mitofinder/meta.yaml)
- Bundled BLAST fails to start on Ubuntu 22.04 (missing `libidn.so.11`); workaround is to point `Mitofinder.config` at system BLAST (https://github.com/RemiAllio/MitoFinder/issues/44)
- Java missing gives "ERROR: java is not installed/loaded" (MiTFi); a failed MiTFi run later shows as `AttributeError: 'bool' object has no attribute 'tRNAs'` in `geneChecker_error.log` (MitoFinder/mitofinder:189-198) (https://github.com/RemiAllio/MitoFinder/issues/14)
- "Gene named COX2 in the reference file(s) are not recognized" appeared with a reference that looked correct (unresolved; only an unrelated `HEG` gene was flagged by the maintainer) (https://github.com/RemiAllio/MitoFinder/issues/40)
- A reference whose feature location joins a part of another accession crashes with "ValueError: Feature references another sequence" (https://github.com/RemiAllio/MitoFinder/issues/66)
- Genes present in a manual BLAST are absent from the result (CO1 missing; BLAST intermediate files empty), maintainer asked for data, no resolution (https://github.com/RemiAllio/MitoFinder/issues/60)
- Circularisation check failed for 30-32 kb genomes with long tandem repeats; the maintainer says the check "is not very sophisticated" (compares only the two ends) (https://github.com/RemiAllio/MitoFinder/issues/62)
- The pillow dependency for the PNG map must be exactly 6.2.2 for Python 2.7 (MitoFinder/README.md:108)
- In the `-a` path with a non-user reference the BLASTX command string lacks a space before `-evalue` (probably dead code because user references are always used) (MitoFinder/geneChecker_fasta.py:134)
- Maintenance: v1.4.2 April 2024, last commit September 2025 (README), 28 open issues, no tagged Python 3 version (https://api.github.com/repos/RemiAllio/MitoFinder)
- Licence: GitHub detects none; the repository has `License/LICENSE` (the bioconda recipe records MIT) and the maintainer says each script carries an MIT text (MitoFinder/License/LICENSE:1) (https://github.com/RemiAllio/MitoFinder/issues/9)

## Limitations by taxon
- Validation in the paper is on ants only; missed genes there were mainly ATP8 (131 of 296 libraries) and the two rRNAs (12S 29, 16S 10) (https://pmc.ncbi.nlm.nih.gov/articles/PMC7497042/)
- Quality depends on how close the reference is; authors write "close reference required" for introns and NUMTs and advise tuning `--nwalk`, E-value, identity and size cut-offs with the reference proximity (MitoFinder/README.md:472-488) (MitoFinder/README.md:534)
- A user with a large termite reference database found it better than MITOS and MitoZ for protein-coding genes and rRNA, but with ARWEN (versions 1.2/1.3) only 16 of 22 tRNAs per genome on average in 40 termite genomes; MiTFi and tRNAscan-SE were added in 1.4 as the answer (https://github.com/RemiAllio/MitoFinder/issues/12)
- tRNA-Trp and tRNA-Cys missing between ND2 and COX1 in a termite compared with MITOS2; maintainer suspected ARWEN (https://github.com/RemiAllio/MitoFinder/issues/8)
- Introns: the default mode does not search for genes with introns (the option `--allow-intron` is required); with it, exons are only similarity merges without boundary search (MitoFinder/README.md:479-482) (MitoFinder/mitofinder:92)
- Octocoral UCE data gave outputs of 15-24 kb with missing genes using a single-species reference (https://github.com/RemiAllio/MitoFinder/issues/51)
- Genomes above 25 kb (default size limit) and genes missing from the reference are not annotated (MitoFinder/mitofinder:869) (MitoFinder/mitofinder:713-716)

## Relevance to our pipeline
- Can MitoFinder annotate a bare assembled mitogenome FASTA with no reads? Yes, in principle: `-a genome.fasta -r reference.gb -o <table> -j ID` skips assembly and goes to contig identification and annotation (MitoFinder/README.md:288); the README test case `Hospitalitermes_medioflavus_NCBI.fasta` is exactly this annotation-only use (MitoFinder/README.md:307). Not run here, so the practical result is untested.
- It cannot work without a GenBank reference of a related mitogenome (`-r`). For a fair benchmark the reference must exclude the genome under test (the README annotation-only test uses the genome's own GenBank record NC_036047 as reference, which leaks the answer) (MitoFinder/test_case/Hospitalitermes_medioflavus_NCBI.fasta:1) (MitoFinder/test_case/Hospitalitermes_medioflavus_NCBI.gb:1) and should be chosen per taxon at the same distance a user would have, e.g. one RefSeq mitogenome of the same family; reference gene names must be standard (COX1...rrnS).
- Settings for a benchmark run: `-o` equal to our `--genetic_code`; `--max-contig-size` above the genome length; check that the input is a single record; `-t mitfi` (default) and also `-t trnascan` or `-t arwen` to see the tRNA effect; compare `*_mtDNA_contig.gb`/`.gff`.
- Expect from design: PCG ends within a 5-codon walk of the BLAST hit, no incomplete-stop (T/TA) convention, no intron boundaries (fragmented cox1/nad5 in cnidarians and annelids with group I/II introns), no rotation to a standard start, tRNAs only as good as MiTFi/ARWEN/tRNAscan-SE. It should therefore be weaker than the MitoZ + `mito_intron.py` + curation chain on intron-containing genomes, and comparable to MitoZ on intron-free ones; this prediction needs the benchmark to confirm.
- Install cost is high: Python 2.7 scripts, system automake/autoconf/gcc/Java, bundled old binaries, bioconda only has 1.4/1.4.1 py27 builds; the realistic route is the Singularity image `library://remiallio/default/mitofinder:v1.4.2` (needs Singularity >= 3, not checked on this machine) or a Python 2.7 conda environment (MitoFinder/README.md:44-51).
- For `--allow-intron` evidence: the existing note in our README that MitoFinder only merges BLASTX pieces and does not search boundaries is confirmed by the MitoFinder README (MitoFinder/README.md:482).
- MitoFinder also uses tRNAscan-SE `-O` with `-X 7` and a genetic-code file; it is a second example of the same mode our `mito_trna.py` runs (see tools/trnascan-se.md).

## References
- Allio R, Schomaker-Bastos A, Romiguier J, Prosdocimi F, Nabholz B, Delsuc F. 2020. MitoFinder: Efficient automated large-scale extraction of mitogenomic data in target enrichment phylogenomics. Mol Ecol Resour 20:892-905 (doi:10.1111/1755-0998.13160)
- Jühling F, Pütz J, Bernt M, et al. 2012. Improved systematic tRNA gene annotation allows new insights into the evolution of mitochondrial tRNA structures and into the mechanisms of mitochondrial genome rearrangements (MiTFi). Nucleic Acids Res 40:2833-2845 (doi:10.1093/nar/gkr1131)
- Laslett D, Canbäck B. 2008. ARWEN: a program to detect tRNA genes in metazoan mitochondrial nucleotide sequences. Bioinformatics 24:172-175 (doi:10.1093/bioinformatics/btm573)
- Code, README and issues (https://github.com/RemiAllio/MitoFinder); singularity container (https://github.com/RemiAllio/MitoFinder_container)
