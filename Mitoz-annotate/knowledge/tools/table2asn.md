# NCBI table2asn (validation of mitogenome submissions)

Every bullet outside "What it does" / "Relevance to our pipeline" ends with a source:
`(path/in/source.py:123)`, `(doi:10.xxxx/...)`, `(https://...)`, or is marked 未确证.
Test inputs and validator output for the experiments below: `knowledge/experiments/table2asn/` (copied 2026-10-20).
Check: `python3 benchmark/check_knowledge.py knowledge/`

## Version checked
- table2asn 1.28.1094 (`table2asn -version`), binary dated 2024-12-14, at `/root/miniconda3/envs/mitointron/bin/table2asn`, checked 2026-10-20 (https://www.ncbi.nlm.nih.gov/genbank/table2asn/)
- On the console the binary prints "This copy of table2asn is more than 1 year old. Please download the current version if it is newer." and still runs (observed 2026-10-20); the MitoZ author states the old `tbl2asn` exits when it prints the same kind of message (https://github.com/linzhi2013/MitoZ/issues/113)
- NCBI states table2asn replaces the "now-obsolete" tbl2asn with very similar operation (https://www.ncbi.nlm.nih.gov/genbank/table2asn/)
- NCBI pages read: table2asn, "Validation Error Explanations for Genomes", "Validation and Discrepancy Report Error Explanations", "Submitting Mitochondrial and Chloroplast Genomes to GenBank" (https://www.ncbi.nlm.nih.gov/genbank/genome_validation/) (https://www.ncbi.nlm.nih.gov/genbank/validation/) (https://www.ncbi.nlm.nih.gov/genbank/organelle_submit/)

## What it does
- Reads a nucleotide FASTA (`.fsa`, with source modifiers in the defline), a 5-column feature table (`.tbl`) and an optional submission template (`.sbt`); builds ASN.1 (`.sqn`), translates every CDS, validates it (`.val`, `.stats`) and can write a GenBank flatfile (`.gbf`). It is the program NCBI runs on submissions, so its messages predict what NCBI curators will bounce.

## Algorithm & key source locations
- Inputs are matched by basename: `.fsa` plus `.tbl` with the same name; `-t` gives the template; `-V v` validates (`.val`, `.stats`), `-V b` writes `.gbf`; `-a s` is "batch set of unrelated sequences"; `-indir`/`-outdir` replaced tbl2asn `-p`/`-r` (https://www.ncbi.nlm.nih.gov/genbank/table2asn/)
- Mitochondrial source modifiers in the defline: `[location=mitochondrion] [topology=circular] [organism=...]` and a genetic code modifier; the NCBI page shows `[gcode=4]`, MitoZ writes `[mgcode=<n>]`, and both the CDS `transl_table` and the record code matter (https://www.ncbi.nlm.nih.gov/genbank/table2asn/) (mitoz/annotate/script/prepare_fsa.py:58)
- Test, `[mgcode=5]` present: a CDS containing an in-frame TGA translates without error; test, `mgcode` absent: same CDS gets `SEQ_FEAT.GenCodeMismatch` "Genetic code conflict between CDS (code 5) and BioSource (code 0)" as a warning, translation still uses the CDS `transl_table` (knowledge/experiments/table2asn/rn/out/nocode.val:7)
- The feature table header `>Feature <id>` (MitoZ spelling) is accepted: the tool logs "Recognized annotation format: five-column feature table" (NCBI documents `>Features <id>`) (https://www.ncbi.nlm.nih.gov/genbank/table2asn/)
- NCBI expectations for metazoan mitochondria: 13 CDS, 22 tRNA, 2 rRNA; CDS rarely overlaps a tRNA by more than a few nt; these CDS "use the mitochondrial TAA stop"; features on both strands except in a few taxa such as mollusks; ND3 of some turtles and birds has a skipped extra base annotated as a join with exception `ribosomal_slippage` (https://www.ncbi.nlm.nih.gov/genbank/organelle_submit/)
- NCBI troubleshooting advice for annotation-program output: partial CDS usually mean the program could not find the endpoints; an essential gene as `misc_feature` or `pseudo` often means an internal stop (or RNA editing); check starts and stops by BLASTP and by alignment with related proteins (https://www.ncbi.nlm.nih.gov/genbank/organelle_submit/)

### Validator messages that matter for mitogenome annotation (exact text from table2asn 1.28.1094 unless noted)
- `SEQ_FEAT.NoStop` ERROR "Missing stop codon": CDS not partial at 3' and not ending in a stop; NCBI fix: extend to the stop, mark 3' partial, or `/pseudo` (https://www.ncbi.nlm.nih.gov/genbank/genome_validation/)
- A CDS ending in `T` or `TA` whose only marker is a `note` ("TAA stop codon is completed by the addition of 3' A residues to the mRNA") still gets `SEQ_FEAT.NoStop`; MitoZ writes exactly this, and its own run on NC_015609.1 gave NoStop for COX2 and ND5 (knowledge/experiments/table2asn/rfix/out/mgcode5.val:7-8) (https://github.com/linzhi2013/MitoZ/wiki/Tutorial)
- Adding `transl_except (pos:<last base(s)>,aa:TERM)` to the same CDS removes the NoStop error: tested for T and TA endings on both strands, e.g. `pos:154..155`, `pos:345`, `pos:complement(382..383)`, `pos:complement(573)`; the gbf then shows the same note and the translation ends at the TERM exception (knowledge/experiments/table2asn/rte/out/x.gbf:32)
- A `transl_except` placed on the wrong bases gives `SEQ_FEAT.CDSrange` (Critical) "Code-break location not in coding region - may be frame problem" plus `SEQ_FEAT.TransLen` ERROR "Coding region extends 2 base(s) past stop codon" (knowledge/experiments/table2asn/rte/out/x.val:1-2)
- `SEQ_FEAT.InternalStop` ERROR "1 internal stops. Genetic code [5]" always comes with `SEQ_INST.StopInProtein` ERROR "[1] termination symbols in protein sequence"; NCBI: check genetic code, CDS location and frame, else `/pseudo` (knowledge/experiments/table2asn/rfix/out/mgcode5.val:5) (https://www.ncbi.nlm.nih.gov/genbank/genome_validation/)
- `SEQ_FEAT.StartCodon` ERROR "Illegal start codon used. Wrong genetic code [5] or protein should be partial" always comes with `SEQ_INST.BadProteinStart` ERROR "gap symbol at start of protein sequence"; tested with code 5: ATG, ATT and GTG starts pass, CTG fails (knowledge/experiments/table2asn/rfix/out/mgcode5.val:6) (https://www.ncbi.nlm.nih.gov/genbank/genome_validation/)
- `SEQ_FEAT.PartialProblemNotSpliceConsensus3Prime` WARNING "3' partial is not at end of sequence, gap, or consensus splice site": produced by MitoZ's `>pos` ends on genes that are not at the sequence end; marking the CDS partial does remove NoStop (knowledge/experiments/table2asn/rfix/out/mgcode5.val:2)
- `SEQ_FEAT.MissingTrnaAA` WARNING "Missing encoded amino acid qualifier in tRNA": a tRNA without a `product tRNA-Xxx`; NCBI says use `tRNA-Xxx` if the amino acid is unknown (https://www.ncbi.nlm.nih.gov/genbank/genome_validation/)
- rRNA without product: this build reports `SEQ_FEAT.rRNADoesNotHaveProduct` WARNING "rRNA has no name"; NCBI's page lists it as `SEQ_FEAT.InvalidQualifierValue: rRNA has no name` and wants names such as "16S ribosomal RNA" (knowledge/experiments/table2asn/rfix/out/mgcode5.val:4) (https://www.ncbi.nlm.nih.gov/genbank/genome_validation/)
- `SEQ_INST.CompleteCircleProblem` WARNING "Circular topology without complete flag set" appears for every circular record and is harmless for CDS checks; MitoZ's own `.val` shows it (https://github.com/linzhi2013/MitoZ/wiki/Tutorial)
- `SEQ_FEAT.UnnecessaryException` WARNING "CDS has exception but passes translation test": added `exception RNA editing` to a CDS that translates cleanly (knowledge/experiments/table2asn/ri2/out/mgcode5.val:7)
- Introns: a 2-exon CDS `join(1..90,124..186)` in a `[location=mitochondrion]` record with non-consensus intron ends (AG...CC) produced no `NotSpliceConsensus` message in this build, so splice-site consensus was not enforced there; one test only, general rule 未确证 (knowledge/experiments/table2asn/ri1/out/mgcode5.val:2)
- `SEQ_FEAT.ShortIntron`: CDS with an intron shorter than 11 bp (usually a frameshift patch); NCBI suggests one gene with `/pseudo`, or `-c s` which keeps the translation with `/artificial_location="low-quality sequence region"` (https://www.ncbi.nlm.nih.gov/genbank/genome_validation/)
- `SEQ_FEAT.PartialProblem` "Start does not include first/last residue of sequence": NCBI says features in sequences with little splicing must be complete or run to the sequence end and be partial; ignorable in eukaryotes at an intron/exon boundary (https://www.ncbi.nlm.nih.gov/genbank/genome_validation/)
- `SEQ_DESCR.OrganismIsUndefinedSpecies` NOTE for placeholder names such as "Test sp." (MitoZ default species name) (knowledge/experiments/table2asn/mitoz_default_species.val:1)

## Defaults that matter
- MitoZ calls `tbl2asn -t template.sbt -a s -p . -V vb`, i.e. validation level "normal" (`-V v`) plus flatfile (`b`), no `-M n` genome tests and no discrepancy report `-Z` (mitoz/annotate/annotation.py:475-477)
- With `-V v` only sequence and feature checks are needed; submitter and publication messages (`GENERIC.MissingPubInfo`, `GENERIC.MissingPubRequirement`) appear when no template is given and are not annotation problems (https://www.ncbi.nlm.nih.gov/genbank/genome_validation/)
- `-c` cleanup letters in this build include `b` basic (default), `e` extended, `f` fix product names, `s` short-intron exception, `x` extend ends to gaps; the NCBI table2asn page documents only `f`, `x`, `d`, `D` (https://www.ncbi.nlm.nih.gov/genbank/table2asn/)
- `-Z` discrepancy report is "recommended only for annotated genome and transcriptome submissions" (https://www.ncbi.nlm.nih.gov/genbank/table2asn/)

## Known bugs / issues
- table2asn does not judge biological correctness: a CDS that starts and ends at valid codons passes even if the boundary is wrong; NCBI relies on curators and tells submitters to check with BLASTP and alignments (https://www.ncbi.nlm.nih.gov/genbank/organelle_submit/)
- The `SEQ_FEAT.NoStop` on T/TA stops is a feature-table convention problem (missing `transl_except`), not a sequence problem; 126 of 128 ERRORs on 31 scale-worm MitoZ records were NoStop (temp_scripts/Mitoz-annotate/README.md:632)
- Old `tbl2asn` as called by MitoZ can stop after one year without update, giving no `.gbf` and a crash in the summary step (https://github.com/linzhi2013/MitoZ/issues/113)
- Message text differs between NCBI documentation and the installed build (rRNA name, `SEQ_FEAT.StartCodon` text); match on codes, not on strings (https://www.ncbi.nlm.nih.gov/genbank/genome_validation/)

## Limitations by taxon
- Cnidaria and Porifera are documented exceptions to the 13/22/2 gene count (https://www.ncbi.nlm.nih.gov/genbank/organelle_submit/)
- Turtles and birds: ND3 extra base needs `join` plus `exception ribosomal_slippage`, otherwise InternalStop or frame errors (https://www.ncbi.nlm.nih.gov/genbank/organelle_submit/)
- Strand: features are expected on both strands "with few exceptions (for example mollusks)", so an all-one-strand annotation outside such taxa prompts a strand check of tRNAs and rRNAs (https://www.ncbi.nlm.nih.gov/genbank/organelle_submit/)
- Mitochondrial group I/II introns inside animal PCGs (anthozoan, some annelid cox1): NCBI's metazoan checklist does not mention them; how curators want them annotated is 未确证 (https://www.ncbi.nlm.nih.gov/genbank/organelle_submit/)
- tRNA introns, trans-splicing and RNA editing: no mention in the animal mitochondrial checklist (the plant and fungal items mention them); handling for animals 未确证 (https://www.ncbi.nlm.nih.gov/genbank/organelle_submit/)

## Relevance to our pipeline
- `ncbi_validate.py` and `mito_curate.py` use `table2asn -V vb` on each record before and after curation; counting by code, not text, is correct because the installed build rephrases some messages.
- `mito_curate.py` P3 adds `/transl_except` plus note to T/TA stops; the tests above show it must be on the last base(s) of the CDS in the right orientation (`pos:complement(...)` on the minus strand) or table2asn raises `CDSrange` and `TransLen` instead of fixing NoStop.
- The defline must carry `[mgcode=<n>]` (taken from `transl_table`), else every CDS adds a `GenCodeMismatch` warning (ncbi_validate.py already does this).
- `PartialProblemNotSpliceConsensus3Prime` after curation means a MitoZ `>`/`<` end remains: the gene needs manual review or a boundary rescue, not suppression.
- Introns added by `mito_intron.py` are not checked for splice consensus by this build in a mitochondrial record, so a clean validation says nothing about intron boundaries; evidence must come from alignment (see README intron validation).
- Placeholder organism ("Test sp.") only produces a NOTE; real species names should still be set before submission.

## References
- NCBI table2asn (https://www.ncbi.nlm.nih.gov/genbank/table2asn/)
- NCBI validation error explanations for genomes (https://www.ncbi.nlm.nih.gov/genbank/genome_validation/)
- NCBI validation and discrepancy report error explanations (https://www.ncbi.nlm.nih.gov/genbank/validation/)
- NCBI Submitting Mitochondrial and Chloroplast Genomes to GenBank (https://www.ncbi.nlm.nih.gov/genbank/organelle_submit/)
- MitoZ wiki Tutorial with real `.val` output (https://github.com/linzhi2013/MitoZ/wiki/Tutorial)
