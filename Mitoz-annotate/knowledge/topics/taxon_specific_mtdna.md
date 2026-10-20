Checked 2026-10-20.

# Taxon-specific features of metazoan mitochondrial genomes that break annotation tools

Scope: what is unusual per phylum, which annotation failures it can cause (tools that assume 13 protein-coding genes, 22 tRNAs, one circular molecule, standard gene order, no introns, no split genes), example RefSeq accessions, and the source.
Conventions: (1) Bullets labelled "Risk:" are inferences from the cited fact for tools such as MitoZ, MITOS2, MiTFi, ARWEN and MFannot; they were not tested here. The cited source supports the fact, not the tool behaviour. (2) Per-record numbers (length, topology, number of CDS/tRNA/rRNA, transl_table, join/intron features) were counted by parsing the GenBank flat file of the named accession from NCBI efetch on 2026-10-20; they describe the submitter/RefSeq annotation, not ground truth. (3) Facts not seen on a fetched page are marked 未确证.

Baseline expected by most pipelines: 37 genes in bilaterians (13 protein-coding, 2 rRNA, 22 tRNA) on one circular molecule (https://pmc.ncbi.nlm.nih.gov/articles/PMC5633667/).

## Cnidaria

- Medusozoa (Scyphozoa, Staurozoa, Cubozoa, Hydrozoa) mtDNA is linear with terminal inverted repeats (TIRs); the review states TIRs are present in all medusozoan mt-genomes (doi:10.1093/gbe/evw195)
- Linear mapping: scyphozoan, staurozoan and hydrozoan mtDNA map to single linear molecules, but Cubozoa have eight mitochondrial chromosomes and Hydra spp. have two near-equal chromosomes; fragmentation arose independently in Cubozoa and Hydridae (doi:10.1093/gbe/evr123)
- Box jellyfish Alatina moseri: 18 unique genes on eight linear chromosomes of 2.9-4.6 kb (29 kb total), with 16 identical 750-nt telomeres forming inverted repeats (doi:10.1093/gbe/evr127)
- A. moseri chromosome content: 13 protein genes, 2 rRNAs, only one tRNA type (trnM-CAU, three identical copies on three different chromosomes), plus polB and orf314 on a chromosome with no other full-length gene (doi:10.1093/gbe/evr127)
- Subtelomeric recombination in A. moseri puts the 5' part of rnl inside the 5' ends of several protein genes and shares the 3' end of cox3 with cob in a different frame, so some annotated protein N-termini are inferred and unverified (doi:10.1093/gbe/evr127)
- Risk: the 3' tail of cox3 matching cob in another frame, and rnl fragments at gene starts, can produce false overlapping or extended CDS calls and spurious rRNA hits inside protein genes (doi:10.1093/gbe/evr127)
- Hydra vulgaris RefSeq is two separate linear records: NC_011221.1 (7,686 bp, 7 CDS, 1 tRNA, 1 rRNA) and NC_011220.1 (8,194 bp, 6 CDS, 2 tRNA, 1 rRNA); neither alone contains the 13 PCGs (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_011221.1)
- Risk: a one-record, one-circular-molecule pipeline reports missing genes (e.g. atp6, atp8, cox2, cox3, nad2, nad5 missing from NC_011221.1) for fragmented genomes, and any check that requires a circular, single-record genome cannot pass (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_011220.1)
- Hydroidolinan hydrozoans have a complete cox1 copy at each end of the linear chromosome, inside the TIRs (doi:10.1093/gbe/evr123)
- Risk: duplicate cox1 at both ends plus an assembler that circularizes or trims terminal repeats gives a doubled or truncated cox1 and an artificial circular contig (doi:10.1093/gbe/evr123)
- Extra ORFs in Medusozoa: polB (B-type DNA polymerase) and a DNA-binding ORF (orf314) occur in cubozoan, scyphozoan, staurozoan and trachyline hydrozoan mtDNA, and are lost in nontrachyline hydroidolinans; the authors propose a linear-plasmid origin (doi:10.1093/gbe/evr123)
- Example: Aurelia aurita NC_008446.1 is linear, 16,937 bp, annotated with 15 CDS (including dnab and orf324), 2 tRNA, 2 rRNA, transl_table 4 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_008446.1)
- Risk: polB/orf314 are non-standard genes absent from default MitoZ/MITOS2 protein sets, so they are left unannotated or mislabelled as unknown ORFs (doi:10.1093/gbe/evr123)
- Medusozoan tRNAs are two genes, trnM(cau) and trnW(uca) (trnW missing in some taxa, e.g. the coronate Linuche unguiculata partial genome and Alatina moseri); cubozoan partial sequences had no tRNA genes in that dataset (doi:10.1093/gbe/evr123)
- Non-myxozoan cnidarians have lost all but one or two mt-tRNA genes; the loss is paralleled by loss of nuclear-encoded mitochondrial aminoacyl-tRNA synthetases (doi:10.1093/gbe/evw195)
- Risk: tRNA finders run on cnidarian mtDNA can only produce false positives beyond tRNA-Met and tRNA-Trp; a benchmark expecting 22 tRNAs counts 20 false negatives that are real absences (doi:10.1093/gbe/evw195)
- Myxozoa (placed inside Cnidaria) retain only eight protein-coding genes plus unidentified ORFs in the reported genomes and have the fastest substitution rates in the phylum (doi:10.1093/gbe/evw195)
- Octocorallia mt-genomes carry mutS (mtMutS, a putative mismatch-repair protein) as a non-standard gene, likely from horizontal transfer from prokaryotes or viruses; example Renilla muelleri NC_018378.1 (circular, 18,643 bp, 14 CDS including muts, 1 tRNA, 2 rRNA) (doi:10.1093/gbe/evw195)
- Risk: mtMutS is not part of the standard 13-PCG reference sets; HMM-based and BLAST-based annotators leave a large unannotated region or call it a hypothetical ORF (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_018378.1)
- Hexacorallia: all hexacoral mt-genomes contain a group I intron in nad5, many have a second intron in cox1, and the nad5 intron encloses between 2 and 15 standard mitochondrial genes (doi:10.1093/gbe/evw195)
- Example: Metridium senile NC_000933.1 (circular, 17,443 bp) has ND5 as join(7027..7743,9425..10510) with the intron 7744..9424 containing ND1 (7967..8950) and ND3 (8957..9313); COX1 is join(3263..4155,5009..5708) with a homing-endonuclease ORF (4310..4981) in its intron (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_000933.1)
- Example: Acropora millepora NC_081453.1 (circular, 18,366 bp) has ND5 as join(2434..3153,15111..16226), so exon 1 and exon 2 are about 12 kb apart with ND1, ND2, ND6, ATP6, ND4, rrn12, COX3, COX2, ND4L and ND3 between them; ATP8 (16359..16577) overlaps the start of COX1 (16559..18160) (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_081453.1)
- Risk: genes nested inside another gene's intron, with the host gene in two distant exons, break pipelines that assume one gene per locus, order-based sanity checks, and per-gene start/stop validation; homing-endonuclease ORFs inside cox1/nad5 introns can be dropped or called as extra CDS (doi:10.1093/gbe/evw195)
- Cnidarian genetic code: the ancestral metazoan code (UGA = Trp) is used by most non-bilaterians, and at least some octocorals no longer use UAG, effectively reverting to the standard code; the RefSeq records above use transl_table 4 (doi:10.1093/gbe/evw195)
- Anthozoan mt-genomes are mapped as monomeric circular molecules, and most cnidarians have relatively low mitochondrial substitution rates (Myxozoa excepted) (doi:10.1093/gbe/evw195)

## Porifera

- Sponge mt-genomes carry atp9 (an ATP synthase subunit absent from bilaterian mtDNA), present in all four classes (Calcarea, Demospongiae, Hexactinellida, Homoscleromorpha) and also in the choanoflagellate Monosiga brevicollis (doi:10.1093/gbe/evw195)
- Example: RefSeq demosponge records annotate atp9 as an extra CDS, e.g. Aplysina fulva NC_010203.1 (circular, 19,620 bp, 14 CDS, 25 tRNA, 2 rRNA) and Geodia neptuni NC_006990.1 (18,020 bp, 14 CDS, 24 tRNA) (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_010203.1)
- Risk: tools whose protein reference set is the 13 bilaterian PCGs will not annotate atp9, and tools with a gene-order or gene-count check will flag a 14th CDS (doi:10.1093/gbe/evw195)
- tatC (a twin-arginine translocase subunit) is found in the homoscleromorph family Oscarellidae; horizontal transfer cannot be ruled out (doi:10.1093/gbe/evw195)
- Introns: demosponges have three different group I introns in cox1 and two group II introns were recently found; some Plakinidae have one or two group I introns in cox1; glass and calcareous sponges have none (doi:10.1093/gbe/evw195)
- tRNA count: most demosponge mt-genomes contain 24-25 tRNA genes, but tRNA loss occurred repeatedly (Keratosa, Plakinidae and possibly some calcaronean sponges); Keratosa keep only tRNA-Trp and tRNA-Met (doi:10.1093/gbe/evw195)
- Example: Plakortis simplex NC_014868.1 (circular, 18,551 bp) is annotated with 14 CDS and only 5 tRNA, versus 25 tRNA in Halichondria panicea NC_040168.1 (19,571 bp) (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_014868.1)
- Calcareous sponges (Calcarea) have multipartite linear mt-genomes: Clathrina clathrus has six linear chromosomes (about 51 kb total), Petrobiona massiliana more than 150 chromosomes (above 77 kb), Leucosolenia complicata at least 50 additional chromosomes (doi:10.1093/gbe/evw195)
- In L. complicata and P. massiliana, each identified coding sequence sits on its own chromosome (doi:10.1093/gbe/evw195)
- Example: Clathrina clathrus RefSeq NC_021113.1-NC_021116.1 are separate linear records of 7.9-9.3 kb each with 1-4 CDS and 0-8 tRNA; NC_021116.1 carries atp9 as its only CDS plus 4 rRNA fragments; Petrobiona massiliana chromosomes are separate GenBank records of about 1 kb (e.g. KT997992-KT997998) (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_021116.1)
- Risk: a per-record pipeline reports almost every gene as missing and every record as incomplete; a pooled run needs the whole set of chromosomes as input (doi:10.1093/gbe/evw195)
- rRNA genes are split across chromosomes in C. clathrus and C. aurea, and calcaronean sponges show widespread U-insertion mRNA editing and tRNA editing (mismatches in the aminoacyl acceptor stem) (doi:10.1093/gbe/evw195)
- Risk: U-insertion editing means the mRNA differs from the DNA, so DNA-based translation shows frameshifts and premature stops; tRNA editing means acceptor-stem mismatches that tRNA finders score as poor (doi:10.1093/gbe/evw195)
- atp8 is lost from calcareous sponge mt-genomes; it was thought absent in glass sponges, although a plausible candidate was recently reported in Oopsacas minuta (doi:10.1093/gbe/evw195)
- Example: Oopsacas minuta NC_027419.1 (circular, 19,042 bp) is annotated with 14 CDS (including atp9 and an atp8 label), 22 tRNA, 2 rRNA and transl_table 5; Aphrocallistes vastus NC_010769.1 (17,427 bp) is annotated with 13 CDS (atp9 present, atp8 absent) and 18 tRNA (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_027419.1)
- Glass sponges use +1 programmed translational frameshifting; in NC_010769.1 the COX3 CDS is join(3097..3512,3514..3880) and ND6 is join(11086..11255,11257..11653), both flagged /ribosomal_slippage (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_010769.1)
- Risk: an annotator that requires a continuous ORF reports a premature stop and truncated cox3/nad6 at the slippage site (doi:10.1093/gbe/evw195)
- Genetic code deviations: AGR = serine in glass sponges, serine or glycine in calcaronean sponges, AUA = methionine in P. massiliana, and UAG = tyrosine in C. clathrus and C. aurea (doi:10.1093/gbe/evw195)
- Risk: with the wrong translation table, AGR/AUA/UAG codons give wrong amino acids, spurious stops (UAG) and failed start-codon checks; these taxa do not fit a single metazoan code setting (doi:10.1093/gbe/evw195)
- Rates: demosponges and homoscleromorphs evolve slowly, glass sponges intermediate, calcareous sponges high to extremely high (doi:10.1093/gbe/evw195)

## Placozoa

- Trichoplax adhaerens mt-genome is 43,079 bp, large for an animal, and holds several large introns with intron-associated ORFs; five substantially different placozoan mt-genomes have been published although only one species is formally described (doi:10.1093/gbe/evw195)
- cox1 is split: seven to eight exons in three segments with different transcriptional polarities, joined by group I intron trans-splicing; a U-to-C edit occurs in cox1 after splicing (doi:10.1093/gbe/evw195)
- Example: Trichoplax adhaerens NC_008151.3 (circular, 43,079 bp) is annotated with 17 CDS, in which cox1 appears only as separate cox1a and cox1b features (no complete cox1), plus orf168a, orf168b, orf175, orf677, 22 tRNA, 3 rRNA, transl_table 4 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_008151.3)
- Risk: a trans-spliced gene on opposite strands cannot be found by any method that predicts one contiguous CDS per gene; cox1 is reported missing and the segments as unknown ORFs (doi:10.1093/gbe/evw195)
- All placozoan mt-genomes contain six to seven group I introns and up to four different group II introns; the large subunit rRNA is encoded in two pieces (doi:10.1093/gbe/evw195)
- Risk: rrnL in two pieces with introns inside genes defeats global covariance-model rRNA annotation and one-gene-one-feature output; intron-encoded ORFs are called as extra CDS (doi:10.1093/gbe/evw195)
- atp8 appears to be lost from placozoan mt-genomes; most placozoan mt-genomes contain 24-25 tRNA genes, more than bilaterians (doi:10.1093/gbe/evw195)
- Risk: a 22-tRNA expectation underestimates tRNA number (duplicates and extra tRNAs get filtered as redundant), and atp8 is reported missing as a true absence (doi:10.1093/gbe/evw195)

## Ctenophora

- Ctenophores have the smallest known animal mt-genomes: Mnemiopsis leidyi 10,326 bp and Pleurobrachia bachei 11,016 bp, with reduced gene content, highly abbreviated rRNAs and very few intergenic nucleotides (doi:10.1093/gbe/evw195)
- All mt-tRNA genes have been lost in Ctenophora; the two tRNA genes reported for P. bachei are annotation artifacts (doi:10.1093/gbe/evw195)
- Example: Mnemiopsis leidyi NC_016117.1 (circular, 10,326 bp) is annotated with 11 CDS, 0 tRNA, 2 rRNA, transl_table 4; P. bachei is GenBank JN392469 (circular, 11,016 bp) (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_016117.1)
- atp6 has been transferred to the nucleus in two ctenophores and atp8 was found in neither the mitochondrial nor the nuclear genome (doi:10.1093/gbe/evw195)
- rRNAs are reduced: in M. leidyi the small and large subunit rRNAs are no more than 368 and 878 bp (doi:10.1093/gbe/evw195)
- Risk: any tRNA finder run on a ctenophore mtDNA produces only false positives (the P. bachei case); covariance-model rRNA annotation misses or mis-bounds the very short rRNAs; atp6 and atp8 are reported missing as real absences (doi:10.1093/gbe/evw195)
- Substitution rates are high to extremely high, which the authors link to selfing hermaphroditic biology and population fluctuations; UAG is reassigned to serine in some ctenophores (doi:10.1093/gbe/evw195)
- Risk: with fast-evolving proteins, homology-based protein finders (GeneWise/HMM/BLAST against bilaterian references) lose sensitivity and boundaries are uncertain; reference selection from within Ctenophora is needed (doi:10.1093/gbe/evw195)

## Platyhelminthes

- Genetic code: NCBI table 9 (Echinoderm and Flatworm; AAA = Asn, AGA/AGG = Ser, UGA = Trp) has systematic range Rhabditophora; table 21 (Trematode; additionally AUA = Met) covers Trematoda; table 14 (Alternative Flatworm; UAA = Tyr) is documented with a comment that evidence for UAA = Tyr in flatworms was not found (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
- Risk: choosing table 5 instead of 9 or 21 mis-translates AAA (Lys in table 5, Asn in 9 and 21) and AUA (Met in tables 5 and 21, Ile in 9), so a reference-protein comparison or stop-codon check differs by code; RefSeq records of trematodes still use transl_table 9 (e.g. Fasciola hepatica NC_002546.1) (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_002546.1)
- Neodermata (parasitic flatworms) RefSeq records are annotated with 12 CDS and no atp8: Schistosoma mansoni NC_002545.1 (circular, 14,415 bp, 23 tRNA), Echinococcus multilocularis NC_000928.2 (13,738 bp, 22 tRNA), Taenia solium NC_004022.1 (13,709 bp), all with transl_table 9 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_002545.1)
- Egger et al. 2017: atp8 had not been confidently annotated in any flatworm; manual annotation found the full metazoan PCG complement including atp8 in Stenostomum sthenum (Catenulida) and Macrostomum lignano (Macrostomorpha), and a later search detected atp8 in all published polyclad mitogenomes where it had not been annotated; the authors conclude atp8 loss within flatworms is restricted to Neodermata (doi:10.1186/s12864-017-3807-2)
- Egger et al. found that manual, but not automated, genome annotation revealed atp8 in the basally branching free-living flatworms (doi:10.1186/s12864-017-3807-2)
- Example: Macrostomum lignano NC_035255.1 (circular, 14,193 bp) is annotated with 13 CDS (atp8 present), 21 tRNA, 2 rRNA, transl_table 9; Stenostomum sthenum NC_035256.1 (16,944 bp) is annotated with 14 CDS (nad3 listed twice, from a duplicated segment containing nad3 and rrnS), 32 tRNA, 3 rRNA, transl_table 5 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_035256.1)
- Contested in Tricladida: the 2026 Polycelis mitogenome paper reports a universal absence of atp8 in four Polycelis species, attributes difficulty to atp8's rapid substitution rate and extreme divergence, and found no transmembrane-domain signal when replicating published atp8-detection methods (doi:10.1186/s12864-025-12467-z)
- Risk: for free-living flatworms atp8 is a small, extremely divergent gene that automated HMM/BLAST pipelines miss; a benchmark that treats absent-atp8 as correct for Platyhelminthes is right only for Neodermata and unsettled for Tricladida (doi:10.1186/s12864-017-3807-2)
- Start codons: in Paratetraonchoides inermis (Monogenea, 14,654 bp, 12 PCGs, 22 tRNA, 2 rRNA, no atp8) cox1 uses GTG and all other PCGs ATG (doi:10.1186/s13071-017-2404-1)
- Genome content: cestode mitogenomes rarely exceed 14,000 bp, M. lignano is the smallest free-living flatworm mitogenome, and the Stenostomum record contains a duplicated segment (nad3 and rrnS) (doi:10.1186/s12864-017-3807-2)
- Example: Schmidtea mediterranea NC_022448.1 (circular, 27,133 bp, 12 CDS, 22 tRNA, 2 rRNA) is about twice the length of Macrostomum lignano NC_035255.1 (14,193 bp) (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_022448.1)
- Tricladida mitogenomes show frequent gene rearrangements, mainly tRNA positional changes, and noncoding regions with tandem repeats and stem-loop structures (doi:10.1186/s12864-025-12467-z)
- Risk: large noncoding/repeat regions and frequent tRNA translocations make gene-order-based or neighbour-based tRNA assignment unreliable in triclads (doi:10.1186/s12864-025-12467-z)

## Nematoda

- Genetic code: NCBI table 5 (Invertebrate Mitochondrial Code; AGA/AGG = Ser, AUA = Met, UGA = Trp) lists Nematoda (Ascaris, Caenorhabditis) in its systematic range; table 14 lists Platyhelminthes and Nematoda in its range (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
- Most nematode mitogenomes contain 12 PCGs because atp8 is missing in all Spirurida and other more derived nematode groups, plus 22 tRNAs and 2 rRNAs (doi:10.1186/s12864-017-4237-x)
- Enoplea retain atp8: Trichinella spiralis mtDNA encodes a putative atp8, absent from all other nematode mtDNAs examined at that time, and its 22 tRNAs include conventional cloverleaf types (doi:10.1093/genetics/157.2.621)
- Example: RefSeq Caenorhabditis elegans NC_001328.1 (circular, 13,794 bp) and Ascaris suum NC_001327.1 have 12 CDS (no atp8), 22 tRNA, 2 rRNA, transl_table 5; Trichinella spiralis NC_002681.1 (16,706 bp) and Trichuris suis NC_017747.1 (14,436 bp) have 13 CDS (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_001328.1)
- Risk: an expectation of 13 PCGs flags most nematode mitogenomes (those of Spirurida and more derived groups) as missing atp8 although the gene is absent, while a 12-gene filter would wrongly drop a real atp8 in Enoplea (doi:10.1093/genetics/157.2.621)
- tRNA arm loss: many nematode mt-tRNAs lack either the D-arm or the T-arm plus variable loop; in Enoplea (Romanomermis culicivorax) the mt genome encodes no canonical tRNA and all 22 are D-, T- or completely armless, with minimal length 45 nt including CCA (doi:10.1093/nar/gky593)
- Most of these armless tRNA sequences are computational predictions; in R. culicivorax several were verified as processed transcripts carrying a post-transcriptionally added CCA end (doi:10.1093/nar/gky593)
- Example: Romanomermis culicivorax NC_008640.1 (circular, 26,194 bp) is annotated with 14 CDS (nad3 listed twice, no atp8), 32 tRNA features and 2 rRNA; Xiphinema pachtaicum NC_033870.1 (12,489 bp) has 18 tRNA and X. rivesi NC_033869.1 21 tRNA (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_008640.1)
- Risk: cloverleaf-trained tRNA finders (MiTFi, ARWEN, ARAGORN) miss or mis-bound D-/T-armless tRNAs, so tRNA recall drops and anticodon calls fail; whether these specific tools recover Enoplea tRNAs is not tested here (未确证)
- Start/stop codons: Globodera ellingtonae uses TTG for cob, cox1, cox2 and nad1, ATT for cox3 and ATA for atp6, and cox1 ends in an incomplete T* stop (doi:10.1186/s12864-016-3047-x)
- Camallanus cotti uses TTG and GTT start codons and a truncated T-- stop on cox3 presumed completed to UAA by 3' A addition to the mRNA (doi:10.1186/s12864-017-4237-x)
- Camallanus cotti also shows annotation ambiguity: GTT start choices change gene length and overlap/IGR calls, and pseudo-gene fragments (pseudo-nad1, pseudo-nad2, pseudo-nad5, pseudo-tRNA) arise from duplicated segments (doi:10.1186/s12864-017-4237-x)
- Risk: ATN-only start-codon rules reject real TTG/GTT/ATT starts, and incomplete T/TA stop codons are seen as missing stops or frameshifts; in T. spiralis all protein genes appear to initiate with ATN and no TTG or GTT starts were found, unlike the TTG/GTT starts inferred for several other nematodes (doi:10.1093/genetics/157.2.621)
- nad4L: a very small gene (234 bp in Camallanus cotti) that overlaps nad6 by 6 bp, one of 13 gene overlaps of 1-11 bp in that genome (doi:10.1186/s12864-017-4237-x)
- Risk: very short ORFs such as nad4L are the first to fall below length or e-value cut-offs, and the 6-bp nad6/nad4L overlap violates no-overlap assumptions (doi:10.1186/s12864-017-4237-x)
- Multipartite mt genomes: Globodera pallida and G. rostochiensis have at least five or six mini-circles of about 6-9 kb, the first circular multipartite metazoan mt genome reported (doi:10.1186/s12864-016-3047-x)
- G. ellingtonae has two circles: mtDNA-I (KU726971, 17,757 bp: cox1-3, nad2-6, nad4l, cob, trnR, trnT) and mtDNA-II (KU726972, 14,365 bp: atp6, nad1, both rRNAs and 20 tRNAs), sharing a ~6.5 kb noncoding region; initial assembly gave two contigs of 12.6 and 13.1 kb (doi:10.1186/s12864-016-3047-x)
- Risk: a single-circle pipeline sees two incomplete genomes, each missing most of the 22 tRNAs or most PCGs, and the shared 6.5 kb noncoding region confuses circularization and coverage-based assembly checks (doi:10.1186/s12864-016-3047-x)
- Bipartite mt genomes also occur in Rhigonematomorpha: Ruizia karukerae has two circular chromosomes (7,659 bp with six PCGs, 2 rRNA, 9 tRNA; and 7,647 bp with seven PCGs and 16 tRNA) that share a 1,010 bp sequence with duplicate cox2 and trnD, trnG, trnH (doi:10.1038/s41598-018-25759-0)
- Enoplea mitogenome architecture is discontinuous: Trichinellidae and Trichuridae keep an almost perfectly conserved architecture, Capillariidae rearrange tRNAs, and Mermithidae and Longidoridae show almost no similarity to the ancestral architecture; skew inversions create strong base-composition bias (doi:10.1186/s12864-022-08607-4)

## Arthropoda

- Sucking lice (Anoplura) have fragmented mt genomes: Pediculus humanus has its 37 genes on 18 minicircular chromosomes of 3-4 kb with one to three genes each; that study found none in the chewing lice and Psocoptera examined and reported recombination between minicircles (doi:10.1101/gr.083188.108)
- The elephant louse Haematomyzus elephantis (Rhynchophthirina) has 33 genes on 10 minichromosomes of 3.5-4.2 kb with 2-6 genes each (doi:10.1038/srep17389)
- The buffalo louse Haematopinus tuberculatus has ten circular minichromosomes of 2.9-5.0 kb carrying 37 genes in total, one to eight genes each plus one large noncoding region (doi:10.3390/ijms232113092)
- Louse and booklouse (Psocodea) fragmentation: about 30 parasitic lice and booklice species have 2-20 minichromosomes, and Liposcelis shows very low conservation of mt gene order relative to other metazoans (doi:10.1186/s12915-021-01218-7)
- Example: Liposcelis entomophila is split in RefSeq into NC_025503.1 (11,599 bp, 11 CDS, 0 tRNA, 0 rRNA, no atp8 or nad4l) and NC_025504.1 (12,675 bp, 1 CDS, 14 tRNA, 2 rRNA), while L. decolor NC_023839.1 is a single circular 14,405 bp record with 13 CDS and 22 tRNA (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_025503.1)
- Risk: a per-record "13 PCGs, 22 tRNAs, 2 rRNAs, circular" completeness check rejects every minichromosome and falsely flags missing genes, while the same species can have fragmented and non-fragmented congeners (doi:10.1186/s12915-021-01218-7)
- Each Haematopinus tuberculatus minichromosome carries one large noncoding region; H. suis and H. apri also have tRNA pseudogenes, and trnP and trnT share identical stretches of 26-27 bp in H. suis, H. apri and H. asini (doi:10.3390/ijms232113092)
- Risk: tRNA finders can report pseudogenes or recombination-shared tRNA copies in louse minichromosome noncoding regions as extra tRNAs (doi:10.3390/ijms232113092)
- Arachnida: T-arms are missing from many of the 22 tRNAs in 4 of 7 arachnid orders, and all opisthothele spiders have some tRNA genes without well-paired aminoacyl acceptor stems; T-arm loss evolved in parallel and once lost is not regained (doi:10.1093/molbev/msn051)
- Mites (Acariformes) show even further tRNA reduction, comparable to Enoplea nematodes (doi:10.1093/nar/gky593)
- Risk: cloverleaf-based tRNA finders (MiTFi, ARWEN, tRNAscan-SE) lose recall on arachnid tRNAs lacking T-arms or acceptor-stem pairing, and a covariance-model score threshold tuned on insects discards them; the effect on specific tools was not measured here (未确证)
- Example: the RefSeq spider Habronattus oregonensis NC_005942.1 (circular, 14,381 bp) is annotated with 13 CDS, 22 tRNA and 2 rRNA, and Limulus polyphemus NC_003057.1 (14,985 bp) with 13 CDS and 22 tRNA, all with transl_table 5 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_005942.1)
- Mite gene order is variable: a 2023 comparison of 23 gamasid mite mitogenomes (10 families) found a novel gene arrangement type in Echinolaelaps echidninus (15,736 bp, 13 PCGs, 22 tRNA, 2 rRNA, control region of 1,561 bp) and different degrees of rearrangement among taxa (doi:10.1017/s0031182023000367)
- Example: Galendromus (Metaseiulus) occidentalis NC_009093.1 (circular, 24,961 bp) is annotated with 21 CDS (10 genes appear twice), 27 tRNA, 2 rRNA and three control-region features; whether this reflects a true duplication or an assembly/annotation artifact was not established here (未确证) (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_009093.1)
- tRNA genes are the most mobile elements in arthropod mt gene order; closely related species differ in tRNA position and tRNAs were excluded from gene-neighbour analyses because of pervasive mobility (doi:10.1186/s12915-021-01218-7)
- Risk: order-based tRNA assignment, or a reference gene order used to resolve ambiguous tRNA hits, mislabels tRNAs in rearranged groups such as lice and mites; a Drosophila-like ground pattern cannot be assumed (doi:10.1186/s12915-021-01218-7)
- Control region: Drosophila melanogaster NC_024511.2 is circular, 19,524 bp, with 13 CDS, 22 tRNA, 2 rRNA and a large annotated A+T region (control region), versus 14,985 bp for Limulus polyphemus NC_003057.1 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_024511.2)
- Bernt et al. note that RefSeq mitogenome annotations, often used as ground truth, are plagued by misassigned reading directions, erroneous gene names and missing or false positive RNA gene annotations (doi:10.1016/j.ympev.2012.08.023)

## Mollusca

- NCBI table 5 (Invertebrate Mitochondrial) lists Mollusca (Bivalvia, Polyplacophora) in its systematic range, and its comments note that several arthropods translate AGG as lysine instead of serine (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
- Doubly uniparental inheritance (DUI) is reported only in bivalves: two mitochondrial lineages coexist, a female-transmitted F-type and a male-transmitted M-type, so one individual (male gonad) yields two divergent complete mt genomes (doi:10.1098/rsob.240186)
- DUI species carry sex-specific ORFs without recognizable homologies, and modified cox2 (extensions, insertions, duplications); the M mtDNA of the venerid T. literatus has an rrnS duplication, and several F or M mtDNAs of the newly sequenced DUI venerids (e.g. both types of T. conspersus and Paratapes undulatus) do not encompass the full 13 PCGs, 2 rRNA and 22 tRNA set (doi:10.1098/rsob.240186)
- Risk: reads from a male can contain two divergent mt genomes, so assemblers may produce chimeras or a second contig; sex-specific ORFs and extended cox2 are non-standard features that become unannotated regions or frameshift/extension calls, and the F-type reference is wrong for the M-type genome (doi:10.1098/rsob.240186)
- In the hermaphroditic freshwater mussel Anodonta cygnea a putative translocated nad5 portion was found, and hermaphroditic-lineage ORFs (H-ORFs) were analysed in several unionid species (doi:10.1186/s12864-018-4583-3)
- Example: Ruditapes philippinarum NC_031332.1 (circular, 22,089 bp) is annotated with COX2 twice, one copy crossing the origin as join(21551..22089,1..1432) and one at 1762..3330, 24 tRNA and no atp8 feature (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_031332.1)
- Genes spanning the circular origin: Magallana (Crassostrea) gigas NC_001276.1 annotates COX3 as join(18084..18224,1..735) over the 18,224 bp origin (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_001276.1)
- Risk: linearized assemblies or annotation tools that treat the sequence as linear split or lose genes crossing position 1 (cox3 here), and a duplicated cox2 or long unassigned regions inflate genome length and gene counts (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_001276.1)
- atp8 was declared missing in several bivalve species, especially in the genus Mytilus, but it was found in Solemya velum, Musculista senhousia, Venerupis philippinarum and Meretrix lamarckii (doi:10.1371/journal.pone.0153631)
- In Meretrix lamarckii atp8 is not less conserved than other mitochondrial genes by Jin-Nei corrected distance, and atp8 and atp6 are not neighbouring in either the F- or M-type mtDNA (doi:10.1371/journal.pone.0153631)
- Missing atp8 is common in bivalves (Mactridae, Arcidae, Pectinidae); Arcuatulinae, Mytilinae and Mytiliseptinae were identified as lacking atp8, and sequences lacking atp8 were annotated manually by scanning intergenic regions; some researchers propose the peptide is dispensable, transferred to the nucleus, or too short and variable to annotate (doi:10.3390/genes14040910)
- Example: Mytilus galloprovincialis NC_006886.2 (circular, 16,744 bp) is annotated with 12 CDS (no atp8), 23 tRNA and 2 rRNA; Placopecten magellanicus NC_007234.1 (32,115 bp) with 12 CDS (no atp8), 32 tRNA and 2 repeat regions; Magallana gigas NC_001276.1 (18,224 bp) with 12 CDS (no atp8), 23 tRNA and 4 rRNA features (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_006886.2)
- Duplicated rRNA: the NC_001276.1 record lists two s-rRNA features (<3786..>4822 and <5820..>7024) and two l-rRNA features (5103..5704 and 17265..17977) with partial (< >) boundaries on both s-rRNA features; whether that is a true rRNA duplication or an annotation artifact was not established here (未确证) (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_001276.1)
- Risk: genomes of 22-32 kb with duplicated tRNAs/rRNAs and extra repeats produce duplicate hits that deduplication logic merges or discards, and a "13 PCGs" expectation reports atp8 missing even when a short divergent atp8 exists (doi:10.1371/journal.pone.0153631)
- Non-bivalve molluscs follow the standard pattern in RefSeq: Katharina tunicata NC_001636.1 (15,532 bp), Octopus vulgaris NC_006353.1, Biomphalaria glabrata NC_005439.1 (13,670 bp) and Haliotis rubra NC_005940.1 each have 13 CDS, 22 tRNA and 2 rRNA with transl_table 5 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_001636.1)

## Annelida

- Annelid mitogenomes in RefSeq use transl_table 5 and the standard gene complement: Lumbricus terrestris NC_001673.1 (circular, 14,998 bp), Platynereis dumerilii NC_000931.1 (15,619 bp), Riftia pachyptila NC_026860.1 (14,987 bp) and Urechis caupo NC_006379.1 (15,113 bp) each have 13 CDS, 22 tRNA and 2 rRNA (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_001673.1)
- Gene order is highly conserved in most annelids (Pleistoannelida: Errantia and Sedentaria), but this ground pattern does not hold in basal-branching lineages: Owenia, Magelona, Eurythoe and Chaetopteridae show substantially different gene orders, and the first annelid mt genomes encoding genes on both strands were reported (doi:10.1016/j.ympev.2015.08.008)
- Syllidae are unusually variable for annelids in gene order, with two trnM copies in Typosyllis antoni and Typosyllis sp. and many noncoding regions in the longest genomes (doi:10.1016/j.gene.2016.08.050)
- Scale worms (Aphroditiformia) mostly have a conserved gene order but deep-sea polynoids show two novel arrangement patterns, which the authors take as evidence against the assumption that Errantia gene order is conserved (doi:10.1016/j.ympev.2018.04.002)
- Siboglinid mt genomes have a conservative gene order with 37 genes (13 PCGs, 22 tRNA, 2 rRNA) (doi:10.3390/genes15010077)
- Risk: a pipeline that uses annelid ground-pattern order to disambiguate tRNAs, to place trnM duplicates, or to call missing genes fails in Syllidae, Chaetopteridae, Oweniidae, Magelonidae and deep-sea polynoids; two trnM copies would be collapsed as redundant (doi:10.1016/j.gene.2016.08.050)
- Reference-record parsing: Owenia fusiformis NC_028712.1 stores gene names only in /note (e.g. note="COX1") and /product with locus_tag AU115_gp01, without a /gene qualifier, so name-matching comparators that read /gene see no genes (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_028712.1)

## Chordata

- Vertebrates use NCBI table 2 (Vertebrate Mitochondrial: AGA and AGG = stop, AUA = Met, UGA = Trp), tunicates table 13 (Ascidian: AGA/AGG = Gly, AUA = Met, UGA = Trp) (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
- Example: Homo sapiens NC_012920.1 (circular, 16,569 bp) has 13 CDS, 22 tRNA, 2 rRNA, transl_table 2; Petromyzon marinus NC_001626.1 and Eptatretus burgeri NC_002807.1 also use transl_table 2; the cephalochordates Branchiostoma lanceolatum NC_001912.1 and B. floridae NC_000834.1 are recorded with transl_table 5, 13 CDS and 22 tRNA (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_012920.1)
- Risk: a single metazoan code setting mis-calls AGA/AGG (stop in vertebrates, Ser in table 5, Gly in tunicates), producing false internal stops or missing stops at gene ends (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
- Control region: the human record annotates a D-loop as complement(join(16024..16569,1..576)), a feature spanning the circular origin, and ATP8 (8366..8572) overlaps ATP6 (8527..9207) by 46 nt and ND4L (10470..10766) overlaps ND4 (10760..12137) by 7 nt (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_012920.1)
- Risk: a control region crossing position 1 is lost or split when the assembly is linearized at a different start, and overlapping gene pairs (atp8/atp6, nad4l/nad4) violate no-overlap assumptions of ORF finders (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_012920.1)
- In mammals tRNA-Ser(AGY) lacks the D-arm, replaced by a single-stranded connector region (doi:10.1093/nar/gky593)
- Risk: cloverleaf tRNA finders miss the D-armless tRNA-Ser(AGY) in mammals unless a model allows arm loss; the effect on specific tools was not measured here (未确证)
- Tunicate (ascidian) mt genomes are hypervariable in almost all genomic features, with extremely high sequence-divergence rates and rampant gene rearrangements even between congeneric and cryptic species; each of six newly assembled ascidian mitogenomes had a novel gene order (doi:10.1093/gbe/evt081)
- Tunicate atp8 was initially considered absent, but the six assembled mitogenomes encode the small atp8 together with the 13 PCGs; an extra trnM gene is present in ascidians (doi:10.1093/gbe/evt081)
- Example: Ciona intestinalis NC_017929.1 (circular, 14,591 bp) is annotated with 13 CDS, 24 tRNA, transl_table 13; Halocynthia roretzi NC_002177.1 (14,771 bp) with 12 CDS (no atp8 feature), 23 tRNA; Botryllus schlosseri NC_021463.1 (14,928 bp) with 13 CDS and 27 tRNA (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_002177.1)
- Risk: reference-based annotation from a vertebrate or amphioxus reference fails to place genes in tunicates (no conserved synteny), the RefSeq Halocynthia record lacks an annotated atp8 although tunicates carry one, and tRNA counts above 22 (duplicate trnM, extra tRNAs) are filtered as redundant (doi:10.1093/gbe/evt081)

## Echinodermata

- NCBI table 9 (Echinoderm and Flatworm Mitochondrial Code: AAA = Asn, AGA/AGG = Ser, UGA = Trp) has systematic range Asterozoa (starfishes), Echinozoa (sea urchins) and Rhabditophora (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
- Example: Strongylocentrotus purpuratus NC_001453.1 (circular, 15,650 bp), Patiria (Asterina) pectinifera NC_001627.1 (16,260 bp) and the crinoid Florometra serratissima NC_001878.1 (16,005 bp) are each annotated with 13 CDS, 22 tRNA, 2 rRNA and transl_table 9 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_001453.1)
- Risk: translating an echinoderm mitogenome with table 5 (AAA = Lys, as in the standard code) instead of 9 (AAA = Asn) changes a common codon, and table 2 would read AGA/AGG as stops instead of Ser (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
- Gene order: the basal gene orders of the three echinoderm lineages (Crinoidea, Ophiuroidea, and the Holothuroidea/Echinoidea/Asteroidea group) derive from a hypothetical ancestral crinoid order by one rearrangement each; genome architecture is highly conserved in Echinoidea, whereas Ophiuroidea show the highest gene-order differences and large amounts of unassigned sequence (UAS) (doi:10.1016/j.ympev.2010.01.035)
- Risk: large unassigned sequences in ophiuroids invite spurious ORF and tRNA hits and a duplication-driven rearrangement pattern confuses order-based gene naming (doi:10.1016/j.ympev.2010.01.035)

## Hemichordata

- Rhabdopleura compacta (Pterobranchia) has a 15,814 bp circular mt genome with typical metazoan size and gene content but the strongest known strand-specific mutational bias among deuterostomes, with a very GT-rich main-coding strand (doi:10.1186/1471-2148-11-134)
- Genetic code: AGG is assigned to Lys in R. compacta, otherwise observed only in some arthropods, and about 20 C-terminal residues are missing from its ATP8 protein (doi:10.1186/1471-2148-11-134)
- NCBI created table 24 (Rhabdopleuridae: AGA = Ser, AGG = Lys, UGA = Trp) for R. compacta and table 33 (Cephalodiscidae: additionally UAA = Tyr) for Cephalodiscidae (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
- Perseke et al. note that the published Saccoglossus kowalevskii sequence (NC_007438) shows a derived gene order and a different assignment of AAA than the echinoderm code (AAA is not Asn) (doi:10.1186/1471-2148-11-134)
- Example: RefSeq Saccoglossus kowalevskii NC_007438.1 (circular, 17,037 bp) and Balanoglossus carnosus NC_001887.1 (15,708 bp) are annotated with 13 CDS, 22 tRNA and transl_table 9 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_007438.1)
- Risk: default echinoderm-code translation of a pterobranch mt genome gives wrong Ser/Lys at AGR and the extreme strand bias and unusual codon frequencies may reduce the sensitivity of homology-based gene calling; whether MitoZ or MITOS2 accept table 24 or 33 is not verified (未确证)

## Chaetognatha

- Chaetognath mt genomes (Spadella cephaloptera, Paraspadella gotoi) were the smallest metazoan mtDNAs known in 2006 and lack atp6 and atp8; a reanalysis found a unique tRNA-Met in both (doi:10.1016/j.mito.2006.07.004)
- Reanalysis with tRNAscan-SE and ARWEN of the five available complete chaetognath mitogenomes predicts additional "cryptic" tRNA genes, but the total never reaches 22 (doi:10.1016/j.compbiolchem.2016.04.007)
- Example: Spadella cephaloptera NC_006386.1 (circular, 11,905 bp) has 11 CDS (no atp6, atp8), 0 tRNA and 2 rRNA annotated; Flaccisagitta enflata NC_013814.1 (12,631 bp) and Paraspadella gotoi NC_006083.1 (11,423 bp) have 11 CDS and 1 tRNA each, all with transl_table 5 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_006386.1)
- Risk: tRNA finders return many extra or overlapping tRNA hits in intergenic regions and 16S rRNA that benchmark references do not contain, so chaetognath false-positive counts are inflated (doi:10.1016/j.compbiolchem.2016.04.007)
- Published chaetognath tRNA counts disagree: earlier work saw only trnMet, a later study found two and four tRNA genes in two mitogenomes, and two conspecific mitogenomes had different numbers (one and two), so the reference tRNA count is unsettled (doi:10.1016/j.compbiolchem.2016.04.007)

## Rotifera and Acanthocephala

- Brachionus plicatilis (monogonont rotifer) has a mitochondrial genome of two circular chromosomes: mtDNA-I (11,153 bp; ATP6, COB, NAD1, NAD2, 13 tRNA and both rRNAs) and mtDNA-II (12,672 bp; COX1-3, NAD3-6, NAD4L and 9 tRNA), with mtDNA-I at four times the copy number and a 4.9 kb region of 99.5% identity shared between them (doi:10.1093/molbev/msn058)
- The two-chromosome structure appears stable in 18 Brachionus strains, and no atp8 is encoded by either chromosome (doi:10.1093/molbev/msn058)
- Example: RefSeq stores them as two records, NC_010472.1 (11,153 bp, 4 CDS, 13 tRNA, 2 rRNA) and NC_010484.1 (12,672 bp, 8 CDS, 9 tRNA, 0 rRNA), both transl_table 5 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_010472.1)
- Risk: each record on its own is missing most PCGs, and the shared 4.9 kb region makes reads map to both chromosomes; completeness filters (13 PCGs, 22 tRNAs, 2 rRNA per record) fail on every bipartite genome (doi:10.1093/molbev/msn058)
- The first complete mitochondrial genome of a bdelloid rotifer, Rotaria rotatoria, is 15,319 bp (doi:10.1186/1471-2164-10-533)
- Acanthocephala lack atp8, commonly have tRNAs with nonstandard structure (lacking the TWC arm or the DHU arm), and Heterosentis pseudobagri has all genes on one strand, several highly divergent hard-to-annotate PCGs, and tRNA genes identified only from a conserved central segment with the anticodon (doi:10.1186/s12864-023-09177-9)
- Example: Leptorhynchoides thecatus NC_006892.1 (circular, 13,888 bp), Paratenuisentis ambiguus NC_019807.1 (13,574 bp) and Macracanthorhynchus hirudinaceus NC_019808.1 (14,282 bp) are each annotated with 12 CDS (no atp8), 2 rRNA and transl_table 5, with 24, 21 and 22 tRNA features respectively (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_006892.1)
- Genetic code: NCBI lists table 14 (alternative flatworm) for Platyhelminthes and Nematoda and table 21 for Trematoda, and the RefSeq Syndermata records above use table 5; no evidence for code 14 or 21 in rotifers or acanthocephalans was found in the fetched sources (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
- Risk: cloverleaf-based tRNA models miss acanthocephalan tRNAs lacking an arm or consisting of the anticodon core only, and automatic annotation also misses divergent PCGs (doi:10.1186/s12864-023-09177-9)

## Xenacoelomorpha

- Xenoturbella (Xenoturbellida) has a conserved gene arrangement: the gene order of X. bocki is most similar to chordates and hemichordates and is inferred as ancestral for deuterostomes (doi:10.1186/1471-2148-9-107)
- Example: Xenoturbella japonica NC_036763.1, X. profunda NC_029220.1, X. monstrosa NC_029219.1 and X. hollandorum NC_029218.1 are circular, 15.2-15.4 kb, with 13 CDS, 22 tRNA, 2 rRNA, transl_table 5 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_036763.1)
- Acoela are different: Symsagittifera roscoffensis lacks atp8, its gene order shows no clear similarity to any other metazoan mitochondrial gene order, cox1 starts with GGT, and trnS1 and trnD lack the DHU-stem (doi:10.1186/1471-2148-10-309)
- Isodiametra pulchra has a larger mitochondrial genome with both rRNAs, 21 tRNAs and only 11 protein-coding genes, a suggested duplicated sequence and many long noncoding sequences, while Paratomella rubra has a compact genome with extensive PCG/tRNA overlap (doi:10.1038/s41598-017-01608-4)
- Example: Symsagittifera roscoffensis NC_014578.1 (circular, 14,803 bp) is annotated with 12 CDS (no atp8) and 21 tRNA; Isodiametra pulchra NC_034948.1 (18,725 bp) with 11 CDS (no atp8, no nad4l) and 22 tRNA (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_034948.1)
- Risk: non-ATN start codons (GGT), missing atp8/nad4l, DHU-stem loss and no reference with similar order make acoel annotation depend on manual curation; the gene-order-based sanity checks cannot be used (doi:10.1186/1471-2148-10-309)

## Bryozoa

- Flustrellidra hispida (Ctenostomata) mt genome is 13,026 bp, one of the smallest metazoan mitochondrial genomes, with a unique gene order; all but one of the 36 genes were identified (trnS2 not found), there are seven tRNA/PCG overlaps (2-11 nt) and seven tRNA/tRNA overlaps (1-8 nt), and nad4, cox2, atp8 and nad3 end in an abbreviated T stop (doi:10.1016/j.ympev.2006.03.007)
- cytb is suggested to terminate on (ACT)AA, which the authors say requires mRNA editing to remove AC; the sequence came from two populations and contains 27 SNP sites in a multi-haplotype contig (doi:10.1016/j.ympev.2006.03.007)
- Example: RefSeq Flustrellidra hispida NC_008192.1 (circular, 13,026 bp) has 13 CDS, 21 tRNA, 2 rRNA; Bugula neritina NC_010197.1 (15,433 bp) and Watersipora subtorquata NC_011820.2 (14,144 bp) have 13 CDS and 22 tRNA, all with transl_table 5 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_008192.1)
- Risk: incomplete T stops inside tightly packed genomes, and genes overlapping tRNAs, lead to missing stop codons, truncated or over-extended gene ends, and tRNAs trimmed or dropped because of the overlaps (doi:10.1016/j.ympev.2006.03.007)

## Brachiopoda

- Lingula anatina (Korea) has a 25,790 bp mt genome composed of 13 PCGs, 2 rRNA and 34 tRNA genes, with dramatic structural changes relative to the previously published Japanese L. anatina record (doi:10.1080/23802359.2017.1407711)
- Example: RefSeq Lingula anatina NC_036679.1 (circular, 25,790 bp) is annotated with 13 CDS, 34 tRNA and 2 rRNA; Terebratalia transversa NC_003086.1 (14,291 bp) and Laqueus rubellus NC_002322.1 (14,017 bp) have 13 CDS and 22 tRNA, all with transl_table 5 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_036679.1)
- Risk: 34 tRNA features and a genome about 1.8 times the usual size violate 22-tRNA and size expectations, duplicated tRNAs may be collapsed as redundant hits, and two intraspecific records disagree in structure so a conspecific reference cannot be assumed correct (doi:10.1080/23802359.2017.1407711)

## Nemertea

- Nemertean mitogenomes are close to the standard pattern: Micrura bella is the largest known for Nemertea at 16,847 bp, encodes the 37 typical metazoan genes with the same arrangement as other Heteronemertea, and has 2,037 noncoding nucleotides at 25 sites including two large noncoding regions of 507 and 508 bp (doi:10.3109/19401736.2015.1060429)
- Example: Lineus viridis NC_012889.1 (circular, 15,388 bp), Cephalothrix simula NC_012821.1 (16,296 bp) and Paranemertes cf. peregrina NC_014865.1 (14,558 bp) are each annotated with 13 CDS, 22 tRNA, 2 rRNA, transl_table 5 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_012889.1)
- Risk: low in the sampled taxa; large noncoding regions of 500 bp are the main item that can attract spurious tRNA or ORF calls (doi:10.3109/19401736.2015.1060429)

## Tardigrada

- Echiniscus testudo (Heterotardigrada) has a 15,817 bp mt genome with 13 PCGs, 20 tRNA, 2 rRNA and an AT-rich control region (78.41% AT), and a mostly conserved PCG order versus the eutardigrade Ramazzottius varieornatus (doi:10.1080/23802359.2018.1495118)
- Example: Hypsibius dujardini NC_014848.1 (circular, 15,169 bp) is annotated with 13 CDS, 23 tRNA and 2 rRNA, and Ramazzottius varieornatus NC_031407.1 (14,428 bp) with 13 CDS, 22 tRNA, both transl_table 5 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_014848.1)
- Risk: tardigrade mt genomes are close to the standard pattern, but the extremely AT-rich genome (78.41% AT in Echiniscus) and the 20-tRNA Echiniscus count mean tRNA recall shortfalls need checking before they are called real losses (doi:10.1080/23802359.2018.1495118)
- Tardigrade-specific genetic code or structural exceptions were not found in the fetched sources (未确证)

## Onychophora

- Epiperipatus biolleyi (Peripatidae) has a gene order considerably different from other arthropods and bilaterians and lacks 9 tRNA genes usually present in bilaterian mt genomes; the missing tRNAs have anticodons for 4-fold degenerate codons while the remaining 13 pair with 2-fold degenerate codons (doi:10.1093/molbev/msm223)
- Opisthopatus cinctipes (Peripatopsidae) in contrast has a highly conserved gene order similar to that of various arthropods, so Onychophora are not uniform (doi:10.1016/j.ympev.2010.05.011)
- Example: RefSeq Epiperipatus biolleyi NC_009082.1 (circular, 14,411 bp) is annotated with 13 CDS, 13 tRNA and 2 rRNA with transl_table 5 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_009082.1)
- Risk: a 22-tRNA benchmark counts 9 false negatives that are real absences, and arthropod-based gene-order or tRNA references misplace genes in Peripatidae (doi:10.1093/molbev/msm223)

## Priapulida

- Priapulus caudatus mt genome is 14,919 bp with few short noncoding regions and several overlapping PCGs and tRNA genes; its gene order is unique among Metazoa but differs from the chelicerate Limulus polyphemus only by a single inversion of 18 contiguous genes bounded by rrnS and trnS(nct) (doi:10.1016/j.gene.2006.10.005)
- Example: RefSeq Priapulus caudatus NC_008557.1 has 13 CDS, 22 tRNA and 2 rRNA, 14,919 bp, and its LOCUS line says linear although the COMMENT states COMPLETENESS: full length; Halicryptus spinulosus NC_020030.1 (14,732 bp) is recorded as circular with the same gene counts (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_008557.1)
- Risk: a topology check against the record flags this complete genome as linear, and compact genomes with PCG/tRNA overlaps trim or drop overlapped tRNAs (doi:10.1016/j.gene.2006.10.005)

## Kinorhyncha

- Echinoderes svetlanae and Pycnophyes kielensis mitogenomes are circular, about 15 kb, with 37 genes all on the major strand, a distinct gene order unique among animals, and duplicated methionine tRNA genes (duplication predates the divergence of the two species) (doi:10.1371/journal.pone.0165072)
- The kinorhynch gene sequences are highly divergent, which precludes accurate phylogenetic inference, and the start codons have a DTD consensus with four types in E. svetlanae and five in P. kielensis (doi:10.1371/journal.pone.0165072)
- Example: Echinoderes svetlanae NC_031873.1 (circular, 15,304 bp) is annotated with 13 CDS, 23 tRNA (22 plus the extra trnM) and 2 rRNA, transl_table 5 (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=NC_031873.1)
- Risk: highly divergent genes weaken reference-based gene calling, non-ATG starts consistent with the DTD consensus are rejected by ATG-only rules, and two trnM copies (23 tRNAs) are collapsed or flagged as duplicates (doi:10.1371/journal.pone.0165072)

## Cross-cutting notes

- RefSeq annotation used as benchmark truth is itself imperfect: it is plagued by misassigned reading directions, wrong gene names, and missing and false positive RNA gene annotations (doi:10.1016/j.ympev.2012.08.023)
- Animal mt-tRNAs have low primary-sequence conservation and "bizarre" structures lacking complete cloverleaf domains, are hard to detect and frequently missed; a re-annotation of 1876 metazoan RefSeq genomes found frequent tRNA gene losses concentrated in basal Metazoa, independent partial losses in Arthropoda, and conserved overlaps of tRNAs in opposite reading directions (doi:10.1093/nar/gkr1131)
- In a comprehensive reannotation of more than 650 bilaterian mitogenomes, Nematoda and Platyhelminthes had extremely long branches and gene rearrangements correlated strongly with long branches (doi:10.1016/j.ympev.2013.05.002)
- Mitochondrial protein-coding gene content shows more variation in non-bilaterian than in bilaterian animals, so a 13-PCG expectation fits bilaterians better than Cnidaria, Porifera, Placozoa and Ctenophora (doi:10.1093/gbe/evw195)
