# tRNAscan-SE 2.x (mitochondrial and organellar modes)

Every bullet outside "What it does" / "Relevance to our pipeline" ends with a source:
`(path/in/source.py:123)`, `(doi:10.xxxx/...)`, `(https://...)`, or is marked 未确证.
Check: `python3 benchmark/check_knowledge.py knowledge/`

Path prefixes used below: `ENV/` = `/root/miniconda3/envs/mitotrna` (the installed copy); `trs/` =
`knowledge/experiments/trnascan` (outputs of the test runs made on 2026-10-20,
copied into the repo; the commands are written out in the bullets so they can be repeated); `Mitoz-annotate/` = our pipeline directory.

## Version checked
- tRNAscan-SE 2.0.13 (bioconda build pl5321hab16a5f_0, Infernal 1.1.5, Perl 5.32.1) installed in env `mitotrna`; GitHub release v2.0.13 published 2026-07-05; checked 2026-10-20 (ENV/conda-meta/trnascan-se-2.0.13-pl5321hab16a5f_0.json:1) (https://github.com/UCSC-LoweLab/tRNAscan-SE/releases)
- Paper: Chan PP, Lin BY, Mak AJ, Lowe TM. 2021, Nucleic Acids Res 49:9077-9096 (doi:10.1093/nar/gkab688)

## What it does
- Finds tRNA genes by covariance-model (CM) search with Infernal and assigns isotype and anticodon. Version 2 added a mitochondrial search for vertebrates (`-M`) and kept a generic organellar mode (`-O`).
- Output: tabular text (Begin/End, type, anticodon, score, note), optional structures, FASTA, BED, GFF3.
- Invertebrate, fungal and protist mitochondria have no dedicated models; `-O` is the only option for them.

## Algorithm & key source locations
- `-M vert|mammal`: Infernal search with one CM per isotype/anticodon (22 vertebrate models plus a Cys model without D-arm), vertebrate mitochondrial genetic code fixed, fast HMM filter (`--mid`) switched on, score cutoff = `organelle_cm_cutoff`; any other model name dies with "Only mammal or vert can be used" (ENV/bin/tRNAscan-SE:1147-1205) (ENV/bin/tRNAscan-SE:1167) (ENV/bin/tRNAscan-SE:1204)
- The branches for invertebrate and fungal mitochondrial models (`-M invert`, `-M fungi`) exist only as commented-out code (ENV/bin/tRNAscan-SE:1192-1202)
- Model files present: `TRNAinf-mito-vert-*.cm` and `TRNAinf-mito-mammal-*.cm` (Ala ... Val, LeuTAA, LeuTAG, SerGCT, SerTGA; vert also Cys-no-darm), no invertebrate set (ENV/bin/tRNAscan-SE.conf:134-135) (ENV/bin/tRNAscan-SE.conf:169-182)
- `-O`: search mode "organelle", Infernal by default (COVE only with `-C`), isotype-specific models disabled, pseudogene filter disabled, score cutoff = `organelle_cm_cutoff` (ENV/bin/tRNAscan-SE:1119-1131)
- In `-O` mode the only CM is the single general model `TRNAinf-1415.cm` (CM name `tRNA1415G`, 90 consensus positions); isotype and anticodon are read from the alignment and the genetic-code file (ENV/lib/tRNAscan-SE/tRNAscanSE/CM.pm:511-516) (ENV/lib/tRNAscan-SE/models/TRNAinf-1415.cm:2)
- `-M` consistency flags written to the Note column: `mito_inconsistent_isotype`, `mito_inconsistent_ac`, `mito_unexpected_ac`, `mito_ac_mislocation` (ENV/lib/tRNAscan-SE/tRNAscanSE/CM.pm:1559) (ENV/lib/tRNAscan-SE/tRNAscanSE/CM.pm:1576) (ENV/lib/tRNAscan-SE/tRNAscanSE/CM.pm:1591) (ENV/lib/tRNAscan-SE/tRNAscanSE/CM.pm:876)
- Isotype naming uses a codon table: the default (cytosolic) table makes an anticodon that pairs with UGA a suppressor ("Sup"); `-g <file>` replaces it with a mitochondrial table (files shipped: vertmito, invmito, ystmito, othmito, echdmito, marsumito, cilnuc) (ENV/lib/tRNAscan-SE/gcode/gcode.invmito:4-6) (ENV/lib/tRNAscan-SE/gcode/gcode.vertmito:4-6) (ENV/bin/tRNAscan-SE:1699-1700)
- Training of the vertebrate models: mitotRNAdb, 1085 vertebrate species grouped into 22 isotype/anticodon sets (mammal set from 282 mammals), two training rounds; a D-arm-less model was built for Ser(GCU); evaluated on 3345 RefSeq vertebrate mitogenomes (doi:10.1093/nar/gkab688)
- Reported vertebrate accuracy: 99.0 % of predictions consistent with RefSeq (72 920), 192 not detected, no false positive at the cutoff in the negative test; ARWEN 87.7 % consistent, MiTFi 98.8 % (doi:10.1093/nar/gkab688)
- Change log items for mitochondria: v2.0.4 `--max` allowed in mito mode; v2.0.7 vertebrate mt-tRNAs lacking the D-arm (other than Ser-GCU) handled and marsupial Asp(GCC) code file; v2.0.10 Cys model without D-arm (armadillos, some reptiles) and `MitoHighConfidenceFilter` for pseudogene/duplicate flags (tRNAscan-SE/change_log.txt:1-30)

## Defaults that matter
- Score cutoff: help text says `-X` default 20, but the config sets 15 for organelle and mito modes (`organelle_cm_cutoff`, `isotype_cm_cutoff.mito_*`); the 2021 paper keeps 20 bits for mitochondrial tRNAs; the measured `-O` default output contains a 16.9-bit hit, so 15 is what runs (ENV/bin/tRNAscan-SE.conf:26-33) (trs/out_O.txt:1) (doi:10.1093/nar/gkab688)
- Lowering the cutoff recovers weak true tRNAs: on Platynereis dumerilii NC_000931.1 (annelid, 15.6 kb, 22 NCBI tRNAs) `tRNAscan-SE -O -X 5 -q --brief` gives 18/22 positions (His 13.7, Leu-TAA 11.1, Ile 10.8, Lys 10.2 bits; Ser-TGA only as Undet/NNN 7.3) against 13/22 at the default, with no extra prediction in this genome (trs/x_O.txt:1) (trs/out_O.txt:1)
- `--max` (no HMM filter) did not add anything for `-O` on that genome (13 hits, same as default) (trs/max.txt:1)
- Speed on the same 15.6 kb genome: `-O` 14.6 s, `-M vert` 1.6 s (measured with `time`); `-M` is fast because it forces the HMM filter (ENV/bin/tRNAscan-SE:1167)
- Coordinates: Begin > End means minus strand; columns are Name, tRNA #, Begin, End, Type, AntiCodon, Intron Begin/End, Inf Score, Note (trs/hdr.txt:1-3)
- Existing output files trigger an interactive prompt; use `-Q/--forceow` in batch jobs (https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/48)
- Config file `tRNAscan-SE.conf` must sit next to the executable (written at install from `.conf.src`); it also fixes `temp_dir: /tmp` (ENV/bin/tRNAscan-SE.conf:1-5) (https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/29)

## Known bugs / issues
- Maintainers state that the current version "was not designed to detect mitochondrial tRNA genes with highly degenerative features such as missing both D-arm and T-arm" and anticodon location may be wrong for them (https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/33)
- `Undet` / `NNN` = the program cannot decide isotype and anticodon because the structure is non-canonical (https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/39)
- `Sup` = suppressor tRNA, anticodon pairs with a stop codon of the codon table in use; with the default (cytosolic) table the mitochondrial Trp(TCA) was called Sup in both genomes tested (code 5 annelid, code 4 cnidarian) unless `-g` is given (https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/35) (trs/out_O.txt:6)
- Measured: `-g ENV/lib/tRNAscan-SE/gcode/gcode.invmito` turns that Sup into `Trp` with the same score 23.1 (trs/og.txt:6)
- Open since 2024: occasional negative coordinates (e.g. `66  -2`) seen in `-G` mode; not reported for `-O`/`-M` (未确证 whether they can occur there) (https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/38)
- Installation: conda `trnascan-se` also installs the internal helper `trnascan-1.4`, which is not the tool to run; compiling from source failed with gcc 15 until v2.0.13 (https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/41) (https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/47)
- Maintenance: active (release 2026-07-05, 10 open issues, GPL-3.0) (https://api.github.com/repos/UCSC-LoweLab/tRNAscan-SE)

## Limitations by taxon
- Mitochondrial models are vertebrate only (mammal subset); maintainers: "-M option currently only supports mitochondrial tRNAs in vertebrates", use `-O` for other organelles (https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/21) (ENV/bin/tRNAscan-SE:1204)
- Paper: version 1.3's organellar mode "was not specifically trained on mitochondrial tRNAs"; its performance was mediocre and it misses mt-tRNAs with highly degenerate structure; future models for insects, nematodes, fungi and other clades are planned because their mt-tRNAs often lack one or two arms (doi:10.1093/nar/gkab688)
- No yeast or protist mitochondrial models; `-O` is the recommended generic fallback, and protist models have low priority because many protists import cytosolic tRNAs (https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/24)
- Plant mitochondria: maintainers also point to `-O` (https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/17)
- Measured, annelid Platynereis dumerilii NC_000931.1 (22 NCBI tRNAs, position match only, isotype not checked): `-O` finds 13, missing Gln, Arg, His, Pro, Ser (both), Leu (one), Ile, Lys; `-M vert` finds 11 (missing Asn, Gln, Arg, Glu, Leu x2, Ser x2, Ala, Ile, Lys); `-M mammal` finds 7 (trs/out_O.txt:1) (trs/out_Mvert.txt:1) (trs/out_Mmammal.txt:1)
- Measured, cnidarian Metridium senile NC_000933.1 (2 NCBI tRNAs, Met and Trp): `-O` finds both (Trp as Sup), `-M vert` finds only Trp (trs/m_O.txt:1) (trs/m_Mvert.txt:1)
- Our batch benchmark with `-O`: recall 61.2 % (cnidaria/sponge/placozoa, 183 tRNAs, 1 false positive) and 46.1 % (RefSeq annelids, 880 tRNAs, 0 false positives), versus 96.7 % and 97.3 % for MITOS2 (Mitoz-annotate/README.md:665-671)

## Relevance to our pipeline
- `mito_trna.py` calls `tRNAscan-SE -O -q --brief -o <file>` without `-g` and without `-X` (Mitoz-annotate/mito_trna.py:134).
- Direct consequence found in this review: its parser skips any row whose type is not one of the 20 amino acids (`c[4][:3] not in AA3`), so every `Sup` and `Undet` row is dropped silently (Mitoz-annotate/mito_trna.py:138-139). In invertebrate and cnidarian mitochondria Trp(TCA) is therefore lost from the tRNAscan-SE method. Passing a genetic-code file (`-g .../gcode.invmito` for table 5, `gcode.othmito` for 4, `gcode.vertmito` for 2, `gcode.echdmito` for 9) restores it; MitoFinder maps tables to these files the same way (MitoFinder/tRNAscanChecker.py:287-310). Not yet changed in our code.
- Recall is low by design for invertebrates, so tRNAscan-SE should stay only one voting family ("independent method family" counts in `mito_trna.py`); it has high precision (0-1 false positives in 1063 tRNAs), so a tRNAscan-SE hit is good supporting evidence, and a miss is not evidence of absence.
- `-X 5` to `-X 10` would raise recall (measured 13 -> 18 of 22 on one annelid with no extra calls), but it was tested on one genome only; the precision at low cutoff on other taxa is untested.
- `-M vert` is not suitable for non-vertebrates (11/22 on an annelid, misses Asn, Glu, Ala, Leu, Ser).
- Use `-Q` if the output path can already exist; our code unlinks the file first (Mitoz-annotate/mito_trna.py:133).

## References
- Chan PP, Lin BY, Mak AJ, Lowe TM. 2021. tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. Nucleic Acids Res 49:9077-9096 (doi:10.1093/nar/gkab688)
- Code and issues (https://github.com/UCSC-LoweLab/tRNAscan-SE), web server and manual (http://trna.ucsc.edu/tRNAscan-SE/)
- Our use (Mitoz-annotate/mito_trna.py:1-40)
