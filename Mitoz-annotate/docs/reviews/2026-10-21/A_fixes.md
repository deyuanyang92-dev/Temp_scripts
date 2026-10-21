# A_fixes - mito_rrna_hmm.py fixes for .review/A_rrna_hmm.md (2026-10-21, TDD, no commit)

Each item: failing test first (checked it failed for the stated reason), minimal fix, pass.

## 1. Robustness (review: "One bad record or a broken model kills the whole run", "search_dir uses glob", "SeqIO.write truncates")
- nhmmer(): no check=True; non-zero exit -> RuntimeError("nhmmer failed (exit N): <stderr tail>"); message names conda env mitotools.
- search_record: nhmmer/model errors per gene -> report row "failed: <msg>" (RuntimeError, CalledProcessError, OSError, ValueError).
- search_dir: rglob("*.gbf") + rglob("*.gb") (04.final_gb/no_cox1_gene/ now searched); per-record and per-file try/except -> "failed:" row; the "file" column is the path relative to the searched directory.
- GenBank written atomically (mkstemp in the same dir, copymode, os.replace; temp removed on error); a write failure turns that file's "added" rows into "failed: not saved (...)".
- Tests: stderr in message, CalledProcessError, failing record does not stop the others, subdirectory + .gb, failure outside search_record, atomic write keeps original, file mode kept.

## 2. Cluster ranking (review: "best_span seeds from the best SINGLE hit")
- clusters(): every unassigned hit (score order) seeds a cluster joined with _joinable; best_span = highest SUMMED score; same-column hits are still never joined.
- best_spans() also returns the runner-up that is a second copy (disjoint in the genome, >= 50 % column overlap with the best, so a fragment over few columns is not a copy); search_record writes the row "second copy at X..Y not added" only when that copy would itself pass the filters; not added as a feature.
- Tests: A1+A2 (130 bits) beat stray B (80 bits) -> 1000..1930; second copy vs fragment; report row.

## 3. Model coverage and partial ends (review: "No model-coverage check")
- LENG read from the .hmm header (cached per file version; missing LENG/file -> "failed:" row).
- Accept when covered columns (hmm_to-hmm_from+1) >= 0.4 x LENG and length >= 200 nt (MIN_RRNA 300 -> 200); "rejected: low coverage" otherwise; report column "coverage" (2 decimals).
- 5' partial when hmm_from > 30, 3' partial when hmm_to < LENG-30: BeforePosition on the genome start / AfterPosition on the genome end for + strand (5' = start, 3' = end); swapped for - strand. Status "added (partial 5')", "added (partial 3')", "added (partial 5' and 3')". gene and rRNA feature carry the same location.
- Tests: 38 % rejected; partial 5'/3'/both on both strands (6 cases); full-length hit no mark; unreadable model.

## 4. Length window (review: "Upper length window = 1.5 x median"; "Hard cap of 2500 nt")
- manifest.tsv columns p5_len, p95_len (5th/95th percentile, linear interpolation); upper limit 1.3 x p95; lower 0.7 x p5 -> "rejected: short". Manifests without these columns keep the old 1.5 x median rule.
- INTERPRETATION (flag for the caller): the lower limit 0.7 x p5 is applied only to spans that are complete at both ends. A partial span is short by definition and is judged by coverage >= 0.4 and >= 200 nt; applying 0.7 x p5 to it as well would make the 0.4 coverage rule (item 3) almost dead (a 45 % hit is always < 0.7 x p5). Pinned by test_partial_span_is_not_held_to_the_complete_gene_lower_limit.
- build --manifest-only: recompute manifest.tsv from the GenBank references, keep .hmm files (n_used = NSEQ from the existing .hmm header; phyla without a .hmm are skipped); no mafft/hmmbuild needed. Reference cap 2500 -> 4000 nt (MAX_REF_LEN). manifest.tsv written via temp + os.replace; read_manifest closes its file.
- Tests: window on both sides, p95 > median case, legacy manifest, manifest-only with a 3500-nt reference kept and a 4500-nt one dropped, percentiles 520/2960.

## 5. Origin-spanning rRNA (review: "Origin-spanning rRNA is mishandled")
- Circular records (annotations topology) are searched as seq + seq[:2000] (ORIGIN_EXT), so a gene across the origin is contiguous for nhmmer.
- best_spans(hits, n): clusters that START at >= n are dropped (they duplicate genes of the first part; also never reported as "second copy"), positions and the copy test use modulo n. Deviation from the wording "drop hits that start at >= L": the drop is applied to CLUSTERS, after joining, so a piece that starts at the origin copy can still join a piece before the origin (your own example: hits L-437..L cols 1-437 and L..L+538 cols 438-975 -> one 975-nt gene; dropping the second HIT first would have kept only 437 nt as a partial). A cluster consisting only of start >= L hits is dropped exactly as specified.
- Span with end > n -> CompoundLocation, + strand [S..n, 1..E]; - strand [1..E, S..n] (transcript order as mito_rescue); partial markers go to the outer ends. Overlap with existing genes modulo n. Report location "4564..5000,1..538+".
- Tests: appended 2000 nt (and none for linear), + / - strand join, pieces + duplicate (one gene, no "second copy"), appended-only duplicates, modulo overlap, partial across origin.

## 6. gene_of aliases (review: "gene_of() and mito_curate recognise different aliases")
- Imports mito_curate.RRNA_NAMES and key (public module-level table, one source of truth for exact aliases: srrna, lrrna, rrn12, rrn16, 12srrna, 16srrna ...); the old free-text substring match (NAMES, small/large) is kept as the second step. Product-only "12S ribosomal RNA" already worked through the substring step (characterisation test, passed before the change); rrn12 / rrn16 / srRNA / lrRNA failed before.
- Coupling note: mito_rrna_hmm now imports mito_curate (and so mito_trna); mito_curate does not import mito_rrna_hmm (no cycle).
- Not changed (out of scope): the substring match is still loose ("rnl" matches gene "trnL" if such a name sat on an rRNA-type feature; gene_of is only called for rRNA features).

## 7. NOTE text (review: "/note ... copied into the submission .tbl")
- Feature note: "annotated by profile HMM similarity (per-phylum rRNA model, {phylum}, {n_used} references)". The "check boundaries manually" hint is the "note" column of the report row (new column, with "coverage"). Report written by write_report (REPORT columns), used by main().
- Test: exact note text, no tool name / instruction in the feature, hint in the row; report file columns.

## 8. Docs (review: docstring numbers and fetch path)
- Module docstring rewritten: search-level 166 -> 175 (rr_detect_summary_dev.txt) is stated as a different clustering without the filters; pipeline result 17 -> 14 missed, fp 9 -> 10, exact 59 -> 59 (pipeline_ab_dev.txt, run before the new rules); "nine divergent phyla 30 -> 38" and "placozoan" dropped; rules updated (best cluster by summed score, second copy report, origin, coverage, partial ends, p5/p95 window, note); usage points to scripts/fetch_rrna_refs.py.
- scripts/fetch_rrna_refs.py (new): argparse -o / --phyla | --all-metazoa (own copy of batch_mitoz.METAZOAN_PHYLA, a test compares the two) / --set TSV (accession + source_accession columns, genus column with --exclude-genera) / --per-phylum 150 / --pool / --retmax / --seed / --pause. Writes <phylum>.gb, refs_<phylum>.txt (accession.version per line) and fetch_log.tsv (query, date, hits, fetched, excluded, kept, seed). Sampling is seeded per phylum (random.Random("seed:phylum")), so it does not reproduce the 2026-10-20 sample of the old script (global seed 7; NCBI content has also changed). Not run against NCBI (as instructed); unit tests use a fake E-utilities transport.
- Tests: load_set, exclusion + seed reproducibility + genus exclusion, main output files, phylum list equals batch_mitoz.METAZOAN_PHYLA.

## 6b. Found while rebuilding: "16S-rRNA" was counted as 12S
- The old free-text matcher used plain substrings: key "s-rrna" (rrnS) matched inside "16S-rRNA", so references named gene "16S-rRNA" went into the rrnS alignment (Rotifera 2, Cnidaria 1 in the shipped references; the rrnS models were built with a 16S sequence in them). The exact-alias step fixes these names; additionally every free-text key must now start at a non-alphanumeric boundary (regex lookbehind), which also stops "rnl" matching inside "trnL". Test: product "mitochondrial 16S-rRNA", gene "trnL", "mt-rns". On all 3134 reference rRNA features gene_of() gives the same result as the alias-only version used to build the library (checked), so the library needs no second build.

## 9. Library rebuilt (build date 2026-10-21, cap 4000 nt)
- `python3 mito_rrna_hmm.py build --gb <scratchpad>/rref/refs/*.gb -o <fresh dir> --cpu 8` (fresh directory, then copied into models/rrna_hmm, so no stale model can survive); 41 models in 21 phyla, same model set as before; manifest has p5_len / p95_len; all 41 rows checked against the .hmm headers (NAME, NSEQ == n_used, ALPH DNA, LENG > 0, p5 <= median <= p95). models/rrna_hmm = 8.0 MB incl. 21 refs_<phylum>.txt accession lists (accession.version of the reference pool of each modelled phylum) and README.md (updated).
- Changed model inputs: Porifera rrnL 80 -> 139 references (median 2389 -> 2483, p95 3006), Cnidaria rrnL 129 -> 138, Placozoa rrnL 19 -> 22, Rotifera rrnL 33 -> 35 / rrnS 34 -> 32 and Cnidaria rrnS 138 -> 137 (relabelled 16S-rRNA).

## Real-tool checks (not unit tests; scratchpad/smoke.py, scratchpad/cli)
Shipped library + real nhmmer on four benchmark genomes with the rRNA features removed (truth = NCBI):
- NC_016117.1 (ctenophore): l-rRNA now 3638..4396 "added (partial 5')" coverage 0.81 (NCBI 3521..4398); before: a 336-nt 38 % hit was added as a complete gene. s-rRNA 3155..3513 (NCBI 3153..3520) is marked partial 5' and 3' only because the Ctenophora rrnS model (p95 572 nt) is longer than the gene (see limits).
- NC_045864.1: s-rRNA 5070..5673 (NCBI 5060..5646, 587 nt) now accepted (was "rejected: long" with 1.5 x median 356); l-rRNA 4161..5048 (NCBI 4067..4999) partial 5'/3'.
- NC_008151.3 (placozoan, split 16S): l-rRNA 26423..29878 partial 5' added, another cluster over the same columns reported as "second copy at 8316..9793" (it is region a of the split gene, not a copy: the second-copy row is only a report line).
- NC_014578.1 (acoel): both hits rejected by coverage (0.10 and 0.37 of the model < 0.4), the 104-nt false rrnS of the old run is gone.
- `mito_rrna_hmm.py search -i dir` (CLI): top-level and no_cox1_gene/ files searched, files rewritten, report with coverage and note columns, partial ends round-trip through SeqIO (`<3155..>3513`).

## Not done / open
- batch_mitoz.py (not touched): `r["status"] == "added"` at ~3164 in run_rrna_hmm (and the failed/added counts near 3141) must become `str(r["status"]).startswith("added")`, otherwise "added (partial ...)" rows are not counted as added; the table there does not write the new "note" column.
- Review findings not in the decision list and therefore not changed: Placozoa/other split-rRNA references (n_multi), stale models on rebuild (a fresh output directory was used instead), phylum/clade plumbing for FASTA input and the silent --core mitoz override (batch_mitoz.py), nhmmer version in the report, join tolerance/column-gap rule, trying the next cluster when the best lies on an existing gene, mito_env.py mitoz_profiles findings.
- Partial flags are heuristic: hmm_from > 30 / hmm_to < LENG - 30 over-call "partial" when the model is longer than the typical gene (Ctenophora rrnS) and the Placozoa rrnL model mixes split-gene regions (p5 410 / p95 3008 nt, so its length window is wide).
- Interpretation made in item 4 (lower length limit applies to complete spans only) and in item 5 (origin duplicates dropped at cluster level, not hit level): see above; each is pinned by a test and easy to revert.
- scripts/fetch_rrna_refs.py was not run against NCBI (fake transport in tests only).
