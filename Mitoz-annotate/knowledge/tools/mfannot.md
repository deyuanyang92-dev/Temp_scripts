# MFannot (BFL-lab)

Every bullet outside "What it does" / "Relevance to our pipeline" ends with a source:
`(path/in/source.py:123)`, `(doi:10.xxxx/...)`, `(https://...)`, or is marked 未确证.
Check: `python3 benchmark/check_knowledge.py knowledge/`

Path prefixes used below: `Mfannot/` = https://github.com/BFL-lab/Mfannot (the single Perl script `mfannot`);
`MFannot_data/` = https://github.com/BFL-lab/MFannot_data; `Mitoz-annotate/` = our pipeline directory.

## Version checked
- MFannot GitHub tag/release 1.37 (latest tag listed), repository HEAD 0990c85 (2023-10-02, "[FIX] tag extraction (#62)"); MFannot_data HEAD c578242 (2023-11-24); checked 2026-10-20 (https://github.com/BFL-lab/Mfannot)
- The script itself reports `Unknown` as version unless run from a git checkout (it calls `git tag | tail -n 1`), and the man page still says "VERSION 1.1" (Mfannot/mfannot:97-105) (Mfannot/mfannot.1:189-190)
- Paper: Lang BF, Beck N, Prince S, Sarrasin M, Rioux P, Burger G. 2023, Front Plant Sci 14:1222186 (doi:10.3389/fpls.2023.1222186)

## What it does
- Perl pipeline (wrapper around Flip, BLAST, HMMER, Exonerate, ERPIN, EMBOSS, Muscle and BFL-lab helper tools) that annotates mitochondrial and plastid genomes from a FASTA or Masterfile; output is a Masterfile plus optional Sequin/tbl.
- Its distinguishing feature is probabilistic intron prediction: group I / group II intron detection and typing (RNAweasel), Exonerate gene models that span introns, intron boundary refinement with type-specific splice rules, mini-exon search, intronic ORF (homing endonuclease / maturase) annotation.
- Built around fungal, protist and plant organelle genomes; bilaterian animals are outside its stated focus.

## Algorithm & key source locations
- Main step order (all in one script): Flip ORFs, BLAST of ORFs vs the protein library, selection of proteins for Exonerate, gene models with introns, gene fusion check, mini-exon HMM search, RNA genes (tRNA, rRNA, rnpB, 5S, tmRNA) by external programs, intron-boundary adjustment, extra genes by HMM, empty/intronic ORFs (Mfannot/mfannot:430-509)
- ORFs: `flip -m -l 40 -g <code>` (minimum 40 aa, first codon of each ORF translated as M if it is a start codon of the chosen code) (Mfannot/mfannot:991)
- Protein search: `blastp -outfmt 5 -matrix BLOSUM62` of the ORFs against the concatenation of all `*.faa` files of the library directory (the man page still says `.pep`); the gene name is parsed from the header `... mt <gene>; ...` (Mfannot/mfannot:1002-1022) (Mfannot/mfannot:600) (Mfannot/mfannot:1180)
- Custom library rule: every FASTA header must contain `mt`, `pt`, `pt_cyano` or `mt_alpha` followed by a gene name and `;` and must be unique, otherwise the script dies (so extra animal references can be added as `>Species mt cox1; ...` in a new `.faa` file) (Mfannot/mfannot:619-620) (Mfannot/mfannot:630-631)
- Hit filter: E-value must be <= 1e-10, and a hit with < 50 % identity in its first HSP must even reach <= 1e-11 (Mfannot/mfannot:1160-1161) (Mfannot/mfannot:284)
- Gene model with introns: Exonerate `protein2genome --minintron 142 --maxintron $MAXLENINTRONS --intronpenalty -35` with a neutral (flat) splice PSSM for both splice sites and `--geneticcode`; the intron minimum 142 is hard-coded in the command, `--minintronsize` is not used there (Mfannot/mfannot:2502)
- Exonerate alignments are kept as exons only when raw score >= 200; weaker alignments become "empty ORFs" instead (Mfannot/mfannot:2839)
- Alignment pieces are chained (SortAndGroupC4Report / DefineBestSol) so one gene model can join exons around a long intron (Mfannot/mfannot:2832-2920)
- Untyped-intron filter: an Exonerate gap that RNAweasel cannot type as an intron is turned back into coding sequence when its length is a multiple of 3 and it has no stop codon (reported as an insertion), and an untyped last intron is always dropped (Mfannot/mfannot:3027-3066)
- Intron length filter: each intron must be >= minintronsize and <= maxintronsize (multiplied by the number of typed introns in a twintron); otherwise the protein is split into separate genes (Mfannot/mfannot:3096-3118)
- Intron typing: RNAweasel through the `RNAfinder` wrapper with ERPIN models; group I uses `Erpin_models/Introns/*.epn` (IA1, IA3, IB, IB3, IC1, IC2, ID, older A/B models), group II uses `domV-mito-521.epn` / `domV-cp-comp.epn`; the paper says RNAweasel uses ERPIN "and more recently Infernal" (MFannot_data/config/.mfannot_external_programs.conf:195-224) (MFannot_data/models/Erpin_models/Introns/domV-mito-521.epn:1) (doi:10.3389/fpls.2023.1222186)
- Boundary adjustment (Adjust_all_intronic_junctions): both ends of an intron are shifted together by k = -6..+6 nt; group I candidates are scored with SpliceI.mod and accepted only if the splice score is >= 19; group II (Adjust_TypeII) rejects shifts that leave an exon <= 5 nt, subtracts |k| from the score, and writes "splice boundaries uncertain" when the best score is < 40 (Mfannot/mfannot:5637-5757) (Mfannot/mfannot:5735) (Mfannot/mfannot:5807-5813) (Mfannot/mfannot:5825-5829)
- Splice-site models: `Splice/SpliceI.mod`, `SpliceII.mod` (type-specific) and `neutral.pssm` (given to Exonerate); the paper explains that Exonerate accepts only one splice matrix, so pooled group I + II patterns give an almost neutral matrix and MFannot re-adjusts afterwards (MFannot_data/README.md:19) (doi:10.3389/fpls.2023.1222186)
- Mini-exons: AnnotateMiniExonsByHMM builds a profile HMM from homologous proteins, enumerates candidate short exons inside an intron (zones at least `minexonsize` = 10 nt long), scores them with hmmsearch, and merges the best one with a flanking exon; the paper says Exonerate has difficulties with exons shorter than about 4 codons (Mfannot/mfannot:3327-3447) (Mfannot/mfannot:4132-4200) (Mfannot/mfannot:4168) (doi:10.3389/fpls.2023.1222186)
- Trans-splicing / split genes: if the Exonerate exons lie on both strands, or an exon falls outside the gene span and a trial protein without it is not better, the gene is not modelled; MFannot only writes the comment "may be either split by linearization of a circular-mapping genome, or it is trans-spliced (annotation in comments)" (Mfannot/mfannot:2578-2583) (Mfannot/mfannot:2686-2697) (Mfannot/mfannot:3159-3178)
- Intronic ORFs: AnnotateIntronicOrfs; endonuclease HMMs (LAGLIDADG, GIY-YIG, HNH, ...) plus `intronic_orfs.pep` (27 sequences); skipped under `--light` (Mfannot/mfannot:6704-6930) (Mfannot/mfannot:7083-7086) (Mfannot/mfannot:186) (MFannot_data/protein_collections/intronic/intronic_orfs.pep:1)
- tRNAs: `RNAfinder -m tRNA` using ERPIN models `Erpin_models/tRNA/{tRNA,tRNAint_cp,tRNAint_mt,yeast}.epn`; the paper says anticodons are assigned with super-wobble rules and that codon reassignment can make tRNA identity wrong (MFannot_data/config/.mfannot_external_programs.conf:149-153) (doi:10.3389/fpls.2023.1222186)
- rRNA: rns and rnl by HMMsearch wrappers (several sub-models, rnl ends only within about +/-50 nt), 5S and rnpB by ERPIN; the paper states rns/rnl are not identified in Bilateria and Euglenozoa (MFannot_data/config/.mfannot_external_programs.conf:75-147) (MFannot_data/config/.mfannot_external_programs.conf:164-185) (doi:10.3389/fpls.2023.1222186)

## Defaults that matter
- Genetic code defaults to 1 (standard); supply `-g` explicitly (4 cnidarians, 5 invertebrates, 2 vertebrates, 9 echinoderm/flatworm); supported ids are those with an `EGC.*` file (0-6, 9-16, 21-23), so tables 24 and 25 do not exist (Mfannot/.mfannotrc:1) (MFannot_data/EGC/EGC.index:1) (Mfannot/mfannot.1:64-81)
- Start codons come from the `Starts =` line of the chosen EGC file (decoded by us): code 4 allows TTA TTG CTG ATT ATC ATA ATG GTG, code 5 allows TTG ATT ATC ATA ATG GTG, code 2 allows ATT ATC ATA ATG GTG, code 1 allows TTG CTG ATG (MFannot_data/EGC/EGC.4:53) (MFannot_data/EGC/EGC.5:65) (MFannot_data/EGC/EGC.2:40) (MFannot_data/EGC/EGC.1:29) (Mfannot/mfannot:347-354)
- `maxintronsize` 3000, `minintronsize` 142, `minexonsize` 10, minimum ORF 40 aa, BLAST E 1e-10, matrix BLOSUM62, empty-ORF minimum 300 nt, overlap cut 30 (Mfannot/.mfannotrc:2-14)
- `maxintronsize` is the Exonerate `--maxintron`; introns longer than this are not found as one intron (the benchmark note below relates to this) (Mfannot/mfannot:667-668)
- Flags: `--partial` skips the mini-exon step, `--light` skips endonuclease and all-gene HMM searches (Mfannot/mfannot:470-475) (Mfannot/mfannot:7086)
- Required environment: `EGC`, `MFANNOT_EXT_CFG_PATH`, `MFANNOT_MOD_PATH`, `ERPIN_MOD_PATH`, `MFANNOT_LIB_PATH`, plus `BLASTMAT` pointing to a directory that contains the PAM70 matrix (the script dies otherwise) (Mfannot/README.md:63-72) (Mfannot/mfannot:272-274)
- Install is Unix only; tested with BLAST 2.2.26/2.2.31+, HMMER 3.2.1, Exonerate 2.2.0, Muscle 3.8.1551, EMBOSS 6.6.0, ERPIN 5.5.4.serv; also needs BioPerl, PirObject, PirModels, Flip, Umac, HMMSearchWC, RNAfinder, Mf2sqn, tbl2asn (Mfannot/README.md:18-56)

## Known bugs / issues
- Docker/FASTA input: several users get "No masterfile found" when giving a FASTA to the Docker image although the man page says FASTA is accepted; open since 2022 (https://github.com/BFL-lab/Mfannot/issues/49) (Mfannot/mfannot.1:158-159)
- Missing tRNAs are the symptom of a broken or unconfigured RNAfinder/ERPIN: with RNAfinder not working the run still finishes "but with no tRNAs"; an older report shows "Use of uninitialized value" floods and no tRNAs (https://github.com/BFL-lab/Mfannot/issues/69) (https://github.com/BFL-lab/Mfannot/issues/6)
- "Cannot read result for 'IntronI'" is the message when the intron-typing commands failed (open, no answer) (https://github.com/BFL-lab/Mfannot/issues/51)
- `Can't use an undefined value as an ARRAY reference` at the external-program stage reported on yeast assemblies and again in 2025, still open (https://github.com/BFL-lab/Mfannot/issues/41) (https://github.com/BFL-lab/Mfannot/issues/69)
- ORFs with non-ATG start codons under non-standard codes were wrong (issue 16, open); a later pull request added ATC as possible start for intronic ORFs (https://github.com/BFL-lab/Mfannot/issues/16) (https://github.com/BFL-lab/Mfannot/issues/55)
- `--T <tmpdir>` breaks step 1 in the Docker image ("Cannot read from masterfile ... Masterfile_copy"), open (https://github.com/BFL-lab/Mfannot/issues/63)
- HMMsearch failure at "Find extra genes with HMM" in the Docker image (https://github.com/BFL-lab/Mfannot/issues/52)
- Licence is stated inconsistently: LICENSE.md is GPL-3.0 but the script header refers to version 2 of the GPL (Mfannot/LICENSE.md:2) (Mfannot/mfannot:32-34)
- Maintenance: last code change 2023-10-02, about 15 open issues/pull requests, 24 stars; RNAweasel repository untouched since 2022-02-05 and warns about hard-coded paths (https://api.github.com/repos/BFL-lab/Mfannot) (https://github.com/BFL-lab/RNAweasel)
- Our benchmark: MFannot missed the large nad5 introns (> 5 kb) in cnidarians (38 of 38 missed); whether that run used the default `maxintronsize` 3000 is 未确证 (Mitoz-annotate/README.md:703)

## Limitations by taxon
- The authors state MFannot is optimised for eukaryotes other than bilaterian animals; they cite MITOS2 for bilaterian animals and say MITOS2 cannot model introns (doi:10.3389/fpls.2023.1222186)
- Protein library: 312 `.faa` gene files; animal sequences are only Homo sapiens (13 genes), Metridium senile (12), a Geodia demosponge (14) and `Trichop.adh.GRS` (12, header abbreviated, presumably Trichoplax) plus the choanoflagellate Monosiga (26); no annelid, mollusc, arthropod or nematode, so similarity-based gene finding for those relies on distant fungal/protist proteins (MFannot_data/protein_collections/identified/cox1.faa:1) (MFannot_data/protein_collections/identified/nad5.faa:1)
- Not all introns are detected and expert curation is required before GenBank submission (Mfannot/README.md:10)
- ERPIN profiles cannot find "derived" introns lacking the typical core elements; rRNA introns are only detected, not placed in a gene model; twintron detection was still planned (doi:10.3389/fpls.2023.1222186)
- Trans-spliced genes and genes split by linearisation are flagged in comments, not modelled (Mfannot/mfannot:2578-2583)
- Mini-exons shorter than about 4 codons are hard for Exonerate; the HMM mini-exon step can still miss them in complex gene structures (doi:10.3389/fpls.2023.1222186)
- Speed: the paper gives 2-10 min for 20-200 kb genomes, versus 1-2 h quoted for a small bilaterian mitogenome with MITOS2 (doi:10.3389/fpls.2023.1222186)

## Relevance to our pipeline
- `mito_intron.py` is a Python port of the intron part of MFannot. Correspondences:
  - library `models/mfannot_proteins` = `MFannot_data/protein_collections/identified`; gene name from the header regex (Mfannot/mfannot:1180)
  - pre-screen (Flip, blastp E 1e-10 / 1e-11, SelectProteinForExo, 142 nt gap rule) = mfannot:962-1480
  - Exonerate command and score >= 200 chaining = mfannot:2502 and 2832-2920; the port uses `--maxintron 5000`, while MFannot's default is 3000
  - RNAweasel/ERPIN typing is replaced by Infernal with Rfam RF00028/RF00029 plus group II domain models and group I subtype models
  - boundary adjustment +/-6, group I >= 19, group II >= 40 = Adjust_all_intronic_junctions (mfannot:5637-6100)
  - mini-exon step is a different implementation (homolog-alignment gap test) of AnnotateMiniExonsByHMM
  - trans-spliced genes are reported, not modelled, same outcome as MFannot
- What matters when debugging: (1) wrong genetic code (MFannot default is 1); (2) intron longer than `maxintronsize`; (3) no homologous protein in the library for the group (annelid, mollusc, arthropod) so the gene may never reach Exonerate; (4) untyped Exonerate gaps are dropped by the in-frame/no-stop rule; (5) tRNAs are never reliable if ERPIN/RNAfinder is not set up. Use real MFannot only as an independent reference for introns, not for tRNAs or rRNAs of bilaterians.
- Local evidence: NCBI animal mitogenomes with introns, 85 records: real MFannot found 85 introns of which 83 have coordinates identical to `mito_intron.py`; on annelids 13 intron coordinates are identical with real MFannot (Mitoz-annotate/README.md:703) (Mitoz-annotate/README.md:708)
- Installing the full MFannot locally is heavy (Perl + PirObject + 7 BFL-lab helper repos + ERPIN); this is why the port exists. The web server https://megasun.bch.umontreal.ca/apps/mfannot/ gives the real output for spot checks.

## References
- Lang BF, Beck N, Prince S, Sarrasin M, Rioux P, Burger G. 2023. Mitochondrial genome annotation with MFannot: a critical analysis of gene identification and gene model prediction. Front Plant Sci 14:1222186 (doi:10.3389/fpls.2023.1222186)
- Prince S, Munoz C, Filion-Bienvenue F, Rioux P, Sarrasin M, Lang BF. 2022. Refining mitochondrial intron classification with ERPIN. Front Microbiol 13:866187 (doi:10.3389/fmicb.2022.866187)
- Lang BF, Laforest MJ, Burger G. 2007. Mitochondrial introns: a critical view. Trends Genet 23:119-125 (doi:10.1016/j.tig.2007.01.006)
- Code (https://github.com/BFL-lab/Mfannot), data (https://github.com/BFL-lab/MFannot_data), RNAweasel (https://github.com/BFL-lab/RNAweasel), RNAfinder (https://github.com/BFL-lab/RNAfinder)
- Our port (Mitoz-annotate/mito_intron.py:1-60)
