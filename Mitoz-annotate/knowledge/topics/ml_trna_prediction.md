Checked 2026-10-20.

# Machine-learning / deep-learning tRNA gene prediction (focus: metazoan mitochondrial tRNAs)

Scope: what exists beyond covariance-model (CM) and heuristic tRNA finders, whether any of it is usable for mt-tRNA annotation in the Mitoz-annotate pipeline, and what training our own model would take. Every bullet carries a source or the marker 未确证. Searches were run on 2026-10-20 in Europe PMC (title/abstract queries combining tRNA with deep learning, neural network, CNN, LSTM, transformer, language model, SVM, HMM, mitochondrial, organelle genome annotation), GitHub repository search, and web search. Search-based negatives ("not found") are 未确证 by nature.

## Classical baselines

- tRNAscan-SE 2.0 (Chan, Lin, Mak, Lowe 2021): CM-based, isotype- and clade-specific models, with a mitochondrial search mode; no machine learning in its methods; GPL v3 (doi:10.1093/nar/gkab688) (https://pmc.ncbi.nlm.nih.gov/articles/PMC8450103/)
- tRNAscan-SE 2.0 mt models: trained on mt-tRNA genes of 1,085 vertebrate species from mitotRNAdb, 22 isotype/anticodon sets, a custom D-arm-less model for tRNA-Ser(GCU); evaluated only on vertebrate RefSeq mitogenomes (https://pmc.ncbi.nlm.nih.gov/articles/PMC8450103/)
- MiTFi (Juhling et al. 2012): 22 hand-curated covariance models (Infernal), one per mt-tRNA family, Leu split into Leu-CUN and Leu-UUR, seed alignments of 33-69 sequences each; re-annotated 1,876 metazoan RefSeq mitogenomes (doi:10.1093/nar/gkr1131)
- MITOS2 (Donath et al. 2019): runs MiTFi for mt-tRNAs plus conflict resolution against protein/rRNA annotations; reference data derived from RefSeq63 (Metazoa) (https://pmc.ncbi.nlm.nih.gov/articles/PMC11571959)
- MITOS2 software licence: MIT (branch mitos2 of the GitLab repo), reference data on Zenodo under CC-BY-4.0 (https://gitlab.com/Bernt/MITOS) (https://zenodo.org/records/2672834)
- ARWEN (Laslett and Canback 2008): heuristic mt-tRNA finder for metazoa, reported to detect close to 100% of previously annotated genes (doi:10.1093/bioinformatics/btm573)
- ARAGORN (2004): heuristic cloverleaf-based tRNA/tmRNA finder, 99% sensitivity on 1,290 eubacterial, eukaryotic and archaeal tRNA genes (doi:10.1093/nar/gkh152)
- Head-to-head on RefSeq (40,521 mt-tRNA annotations, 1,876 genomes): sensitivity/precision tRNAscan-SE 0.898/0.981, ARWEN 0.977/0.869, MiTFi 0.986/0.979 (https://pmc.ncbi.nlm.nih.gov/articles/PMC3326299)
- Head-to-head on 3,345 vertebrate mitogenomes (73,674 annotated mt-tRNAs): agreement with RefSeq tRNAscan-SE 2.0 99.0%, MiTFi 98.8%, ARWEN 87.7%; false positives in 17,500 virtual genomes: tRNAscan-SE 2.0 none, MiTFi none below its recommended e-value, ARWEN 34,634 (https://pmc.ncbi.nlm.nih.gov/articles/PMC8450103/)

## ML / DL methods

Summary: only one deep-learning method written specifically for tRNA gene sequences was found (tRNA-DL), plus one older ensemble classifier; neither is mitochondrial-specific, neither is a genome scanner, and no published method comparing ML against MiTFi/MITOS2 on metazoan mt-tRNAs was found (未确证, search-based; Europe PMC REST queries run 2026-10-20, https://www.ebi.ac.uk/europepmc/webservices/rest/search).

### tRNA-DL (Gao, Wei, Hakonarson; Hum Hered 83(3):163-172)

- Year/venue: Europe PMC lists 2018, the NJIT repository lists 2019 (doi:10.1159/000493215)
- Task: binary classification of a candidate sequence as true tRNA versus tRNAscan-SE false positive found in coding sequences; it is a post-filter to be coupled with tRNAscan-SE, not a stand-alone scanner (doi:10.1159/000493215)
- Input: one-hot encoded nucleotide sequence, no handcrafted features; the repo's selected best model takes input_length=134 with 4 channels (https://github.com/stella-gao/trna-dl)
- Model: 13 architectures compared (Conv1D, Conv2D, LSTM, BiLSTM, BiGRU combinations); selected best "CCMBLF" = 2x Conv1D + MaxPool + bidirectional LSTM + dense layers (https://github.com/stella-gao/trna-dl)
- Training data: positives and negatives in the repo fasta/ folder; I counted 426 positive sequences (headers tdbR..., apparently tRNAdb IDs) and 1,183 negative sequences (ENA "PREDICTED PROTEIN" records) after CD-HIT (https://github.com/stella-gao/trna-dl/tree/master/fasta)
- Mitochondrial tRNAs in training data: 未确证 (the paper is paywalled; the abstract and repo do not say, and positives are only labelled by tRNAdb-style IDs)
- Reported accuracy: the abstract states the best model outperformed existing methods and reduces tRNAscan-SE false positives; no numbers were retrievable and MiTFi/MITOS2/ARWEN are not named in the abstract (doi:10.1159/000493215)
- Accuracy versus MiTFi/MITOS2 on mt-tRNAs: 未确证 (no evidence found that it was ever tested)
- Code: GitHub repo with scripts and trained .hdf5 models, no licence file (GitHub API reports none), last push 2018-06-15, Keras-1-style code (border_mode, nb_filter arguments); repo identity as the authors' official release is inferred from the README repeating the abstract (https://github.com/stella-gao/trna-dl)
- Runnable offline: models are shipped in model/, so inference needs no network once a legacy Keras environment exists; not tested by me, so 未确证
- Maintained: no commits since 2018-06-15 (https://github.com/stella-gao/trna-dl)

### tRNA-Predict (Zou et al. 2015, Mol Inform)

- Method: ensemble classifier "LibMutil" over three handcrafted feature sets, to separate functional tRNAs from tRNAscan-SE false positives (doi:10.1002/minf.201500031)
- Training data: 623 tRNA sequences from tRNAdb 2009 as positives, tRNAscan-SE false positives as negatives; 95.1% accuracy by 10-fold cross-validation (doi:10.1002/minf.201500031)
- Mitochondrial coverage and comparison with MiTFi: 未确证 (not stated in the abstract)
- Availability: web server only (datamining.xmu.edu.cn/~gjs/tRNA-Predict); my request on 2026-10-20 returned HTTP 400, no source code found, licence unknown 未确证 (http://datamining.xmu.edu.cn/~gjs/tRNA-Predict)

### Historical: neural network on tRNA gene sequences (Sun et al. 1995)

- Used an artificial neural network to quantify similarity among tRNA gene sequences and to flag a new sequence as tRNA-like; not a genome scanner, no code (doi:10.1089/cmb.1995.2.409)

### General ncRNA / RNA language models (not tRNA gene finders)

- RNA-FM (arXiv 2022): 12-layer BERT encoder, masked-token pretraining on 23.7 M RNAcentral100 ncRNA sequences; downstream tasks are structure and function prediction, not tRNA gene detection; MIT licence, weights on Hugging Face, last push 2025-05-27 (https://github.com/ml4bio/RNA-FM) (https://arxiv.org/abs/2204.00300)
- RiNALMo (Nat Commun 2025): 650 M-parameter RNA language model pretrained on 36 M ncRNA sequences; the repo ships a fine-tuned head for secondary structure of ArchiveII tRNA sequences (structure of a given tRNA, not gene detection); Apache-2.0 code, weights on Zenodo, last push 2026-10-08 (https://github.com/lbcb-sci/RiNALMo) (doi:10.1038/s41467-025-60872-5)
- Whether either model has been tested for locating tRNA genes in a genome or for mitochondrial tRNAs: 未确证 (not found in the abstracts or READMEs I read)
- ncRNA family classifiers benchmarked by Creux et al. 2024 (nRC, RNAGCN, ncrna-deep, MFPred, ncRDense, NCYPred): they classify given sequences into classes and the paper states genome scanning is not addressed; whether tRNA is among the 13 classes was not confirmed from the text I read 未确证 (doi:10.1371/journal.pcbi.1012446)
- MuStARD (Georgakilas et al. 2020): multi-branch CNN that scans genomic regions for small-RNA loci, demonstrated on pre-miRNA and snoRNA, not on tRNA; code on GitLab, last activity 2020-03-24 (doi:10.1038/s41598-020-66454-3) (https://gitlab.com/RBP_Bioinformatics/mustard)
- Neural-SCFG (student project, 2026): transformer encoder + stochastic context-free grammar for tRNA secondary structure; 4,054 GtRNAdb sequences from 16 bacterial and mammalian genomes; median F1 0.358; not peer reviewed, no licence, not mitochondrial, not gene detection (https://github.com/johnmccausland/Neural-SCFG)
- PON-mt-tRNA (2016): machine-learning classifier of human mt-tRNA variant pathogenicity, not gene detection (doi:10.1093/nar/gkw046)
- Names suggested in the brief, tRNAdetector and DeepTRNA: no tool of these names found in Europe PMC, GitHub or web search on 2026-10-20 (未确证)
- Organelle genome annotation by deep learning (title/abstract search for mitogenome/organelle annotation with deep learning, neural network, language model): 0 hits in Europe PMC on 2026-10-20 (未确证)

## Training data available

- RefSeq metazoan mitogenomes: 1,876 genomes / 40,521 mt-tRNA annotations (RefSeq v39, Juhling 2012); 3,345 vertebrate genomes / 73,674 mt-tRNAs (tRNAscan-SE 2.0 evaluation); 13,387 metazoan mitogenomes as of Sep 2022 (https://pmc.ncbi.nlm.nih.gov/articles/PMC3326299) (https://pmc.ncbi.nlm.nih.gov/articles/PMC8450103/)
- Licence of RefSeq data: NCBI states it places no restrictions on use or distribution, while submitters may claim rights (https://www.ncbi.nlm.nih.gov/home/about/policies/)
- Label quality warning: Juhling 2012 disregarded the Ser and Leu distinctions in evaluation because RefSeq shows many misannotations of that type (https://pmc.ncbi.nlm.nih.gov/articles/PMC3326299)
- MITOS2 re-annotation of all 13,387 RefSeq metazoan mitogenomes: 287,415 mt-tRNA genes, about 95% (274,975) cloverleaf; raw results (zip files, 18 MB) on Zenodo under CC-BY-4.0; labels are model-derived, not curated (https://zenodo.org/records/13789952) (doi:10.1093/gbe/evae232)
- The Ozerova text counts 13,387 mitogenomes and 287,415 tRNAs while its abstract says more than 12,500 and more than 250,000; use the main-text numbers with care (doi:10.1093/gbe/evae232)
- Hand-curated aberrant sets: custom armless CMs and nematode literature-based CMs on the same Zenodo record; hand-curated Mermithidae alignments from Juhling 2012b retrieve an almost complete tRNA set in Romanomermis (https://zenodo.org/records/13789952) (doi:10.1093/gbe/evae232)
- MITOS2 reference data (Zenodo concept DOI 10.5281/zenodo.2672834, CC-BY-4.0): refseq63m.tar.bz2 11.3 MB and refseq89m.tar.bz2 19.6 MB include automatically derived covariance models and fasta databases (https://zenodo.org/records/2672834)
- MiTFi seed alignments and models: distributed with MiTFi (also in the MITOS repo, directory mitos/mitfi); the project URL given in the paper, www.bioinf.uni-leipzig.de/software.html, returned 404 on 2026-10-20; separate licence for the MiTFi models 未确证 (https://gitlab.com/Bernt/MITOS)
- mitotRNAdb: 1,085 vertebrate species were used by tRNAscan-SE 2.0 for its mt models; the server http://mitotrnadb.bioinf.uni-leipzig.de/ returned HTTP 502 on 2026-10-20; licence and invertebrate coverage 未确证 (https://pmc.ncbi.nlm.nih.gov/articles/PMC8450103/)
- tRNAdb 2009: more than 12,000 tRNA genes classified by amino acid specificity, NCBI taxonomy integrated; the host trna.bioinf.uni-leipzig.de answered HTTP 200 on 2026-10-20; licence 未确证 (doi:10.1093/nar/gkn772)
- GtRNAdb (nuclear): Release 22, September 2024, 4,867 genomes (4,038 bacterial, 609 eukaryotic, 220 archaeal); the front page states no licence and does not say whether mt tRNAs are included, so mitochondrial content and terms are 未确证 (https://gtrnadb.ucsc.edu/)
- GtRNAdb 2.0 description: more than 367,000 predictions in more than 4,370 genomes at publication (doi:10.1093/nar/gkv1309)
- tRNA-DL benchmark sets (426 positives, 1,183 negatives after CD-HIT, tRNAdb IDs / ENA coding-sequence false positives) are public in the repo but small and not mitochondrial-specific as far as I could tell (https://github.com/stella-gao/trna-dl/tree/master/fasta)

## Assessment for our pipeline

- Is there a published ML method that beats MITOS2/MiTFi on metazoan mt-tRNAs? None found; the only tRNA-specific DL paper (tRNA-DL) targets tRNAscan-SE false positives, reports no numbers in its abstract, and shows no comparison with MiTFi or ARWEN (未确证, absence of evidence from the searches above; doi:10.1159/000493215)
- Headroom on canonical vertebrate-like mitogenomes looks small: MiTFi precision 0.979 and tRNAscan-SE 2.0 with zero false positives in 17,500 simulated genomes leave little to gain from a filter, caveat that RefSeq labels were partly produced by similar tools (circular agreement) 推测 / 未确证 (https://pmc.ncbi.nlm.nih.gov/articles/PMC3326299) (https://pmc.ncbi.nlm.nih.gov/articles/PMC8450103/)
- Where the real gap is: armless and truncated mt-tRNAs. MITOS2 uses only cloverleaf models, so aberrant tRNAs are missed or over-predicted with a false extra arm; truncated custom CMs raise sensitivity but D-arm loss is still not detected and manual checks are needed; clade-specific hand-curated CMs recover almost everything in Romanomermis (doi:10.1093/gbe/evae232)
- Aberrant mt-tRNAs are widespread, not rare: Acanthocephala, Nematoda, Acariformes and Araneae are hotspots, mt-tRNA-Ser is most often affected (4,339 mitogenomes) (doi:10.1093/gbe/evae232)
- Judgement: a DL model would have to learn armless shapes from very few curated examples (a handful of clades), so generalisation to unseen clades is unlikely to beat clade-specific CMs; the cheaper route is adding clade-specific or truncated CMs and a rescoring step 推测 / 未确证
- Requirement 1, data: positives from the MITOS2 re-annotation (287,415 genes) carry the biases and errors of the CM method that produced them, so a model trained on them would mostly reproduce MiTFi; a trustworthy test set needs independent curation (experimentally supported or manually curated genomes) 推测 / 未确证
- Requirement 2, leakage control: split by genus or, better, family/order, with near-identical sequences removed, because mt-tRNAs within a genus are near-identical (MiTFi seed construction already excluded nearly identical sequences, 33-69 sequences per family) 推测 / 未确证 (https://pmc.ncbi.nlm.nih.gov/articles/PMC3326299)
- Requirement 3, negatives: windows from the same mitogenomes outside annotated tRNAs (rRNA, protein-coding, control region, both strands), dinucleotide-shuffled sequences, and hard negatives from lax-threshold ARWEN/MiTFi false positives (ARWEN produced 34,634 on virtual genomes); tRNA-DL used tRNAscan-SE false positives in coding sequences as negatives 推测 / 未确证 (https://pmc.ncbi.nlm.nih.gov/articles/PMC8450103/) (doi:10.1159/000493215)
- Requirement 4, armless tRNAs: need explicit labels for four shapes (cloverleaf, D-arm loss, T-arm loss, armless) and a boundary tolerance, since the 3' CCA is added post-transcriptionally and boundaries are model-defined; the shapes are only curated for a few clades (doi:10.1093/gbe/evae232)
- Requirement 5, anticodon and isotype labels: trnL1/trnL2 and trnS1/trnS2 cannot be assigned from the anticodon alone in the MiTFi work (Leu first modelled as one set, then split Leu-CUN/Leu-UUR; Ser groups UCN vs AGY), and RefSeq has many Ser/Leu misannotations, so these labels are noisy; gene-order or clade context features would be needed 推测 / 未确证 (https://pmc.ncbi.nlm.nih.gov/articles/PMC3326299)
- Likely gain from training our own: small for canonical tRNAs, possibly useful as a false-positive rescoring layer or for isotype disambiguation, unproven for armless tRNAs; I would not expect a publishable improvement without a new curated benchmark 推测 / 未确证
- Reusing existing weights is not a short cut: tRNA-DL is a 2018 Keras-1-era filter without a licence, general RNA language models (RNA-FM MIT, RiNALMo Apache-2.0) would need task-specific fine-tuning and curated mt data first, and were not evaluated on tRNA gene finding 推测 / 未确证 (https://github.com/stella-gao/trna-dl) (https://github.com/ml4bio/RNA-FM) (https://github.com/lbcb-sci/RiNALMo)
- Suggested next step: keep MITOS2/MiTFi as the primary tRNA caller and add a held-out, genus-split benchmark (including armless clades) before testing any ML rescoring; judge ML only by head-to-head recall/precision on that set 推测 / 未确证
