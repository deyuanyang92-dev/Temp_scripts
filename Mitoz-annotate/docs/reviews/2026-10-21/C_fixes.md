# Review C fixes (mito_curate.py, ncbi_validate.py; TDD, one test per finding)

## Finding 1 (Critical) ncbi_validate.feature_lines: partial markers on minus strand
- Tests (tests/test_ncbi_validate.py): `test_tbl_partial_markers_first_column_lt_second_gt_on_minus_strand`
  (FAIL before: `>150\t40`), `test_tbl_partial_markers_plus_strand_unchanged`,
  `test_minus_strand_partial_cds_survives_table2asn` (FAIL before: CDS dropped from rec.gbf; skipped without table2asn).
- Fix: strand -1 parts: first column = end with `<` if AfterPosition(end); second column = start with `>` if BeforePosition(start).

## Finding 15 strip internal notes from the .tbl
- Test `test_tbl_skips_internal_review_notes` (FAIL before: both tool notes written).
- Fix: `INTERNAL_NOTE = r"check .* manually|mito_\w+\.py"`; matching `note` values are skipped in feature_lines only
  (the .gbf keeps them). Other notes (e.g. "start codon not determined") stay.

## Finding 6 per-part neighbour spans
- Test `test_wrapped_trna_does_not_block_p1` (FAIL before: no P1 row, start stayed at 307; control without the tRNA moved to 301).
- Fix: new `spans(f)` (1-based per part). `upstream_end`, `downstream_trna_start`, rrna() neighbours and the rRNA/CDS
  overlap (R3) loop over parts. In `upstream_end` the "different frame" skip applies only to single-part CDS (a join part's
  own start does not give its frame; joins therefore always block, conservative).
- `pcg_overlaps` only ever compares single-interval (`simple`) CDS, so there is no hull there; left as is.

## Finding 7 R1 inside gene
- Tests `test_rrna_not_trimmed_for_inside_gene_when_loss_large` (FAIL before: rRNA 1471..2300), `test_rrna_trimmed_for_inside_gene_near_the_end` (R1 still trims, 1001..2200).
- Fix: loss = nt discarded with the gene (gene start to rRNA end, or rRNA start to gene end, the shorter side kept as before).
  Trim only if loss <= 10 % of the current rRNA length; else keep the rRNA and log one R3 row per gene
  "gene X inside rRNA, not trimmed (would lose N nt) - check manually (spurious tRNA call?)". Existing
  `test_rrna_never_contains_a_trna` (loss 8 %) stays green.

## Finding 9 (+C4) never leave CDS without /translation
- Tests: `test_translation_kept_for_every_cds` (KY753833 through main(): 13 in, 0 out before; 13/13 after, all start with M, equal to extract().translate()),
  `test_translation_alternative_start_is_methionine` (ATT start, + / - / join, table 5 -> `M...`),
  `test_translation_not_written_with_internal_stop_and_partial_start_plain` (internal stop -> none; 5'-partial -> plain `I...`).
- Fix: new `Curator.translation(f)` called for every CDS at the end of run(): `f.extract(rec.seq)`, drop codon_start-1,
  trim to a multiple of 3, translate with the record's genetic code, drop one trailing `*`, first codon in the table's start
  codons and not 5'-partial -> `M`; internal `*` -> /translation not written. The older `pop("translation")` calls in pcg() /
  across_origin() stay (they only guard an aborted run).

## Finding 10 P6 stop test in frame only
- Test `test_p6_stop_tested_in_frame_only` (ATG + GCC*40 + ATAA, 127 nt: no P6 before; P6 + P4 after).
- Fix (`pcg_check`): `has_stop = last 3 nt in stops` only if len % 3 == 0, else tail in (T, TA). A tail without T/TA also logs
  P4 "length N is not a multiple of 3 and there is no T/TA tail ...". P6 uses `has_stop`.

## Finding 11 O1 cap
- Test `test_origin_extension_capped_by_feature_at_position_1` (FAIL before: O1 extended 152 nt into the feature at position 1).
  Cases: stop 17 nt in -> O1; stop 152 nt in -> no O1, P6 partial; no feature at position 1 -> 300 nt limit still applies.
- Fix: `MAX_ORIGIN_OVERLAP = 60`. The continuation starts at genome position 1 (+ strand) or L (- strand, transcript wraps to
  the genome end); if another CDS/tRNA/rRNA (any strand, per part) covers that position the search window is 60 nt, else 300.
  (Generalised from "position 1 of the record": for a minus-strand CDS the continuation lies at the other end of the record.)
  Not done: start-codon / internal-stop checks on the joined CDS (outside decision list).

## Finding 12 PCG loop before R1
- Test `test_cds_extension_does_not_overlap_following_rrna` (FAIL before: CDS extended to 1531, rRNA still starting 1454).
- Fix: run() order is now names -> CDS loop (pcg, across_origin, pcg_check) -> pcg_overlaps -> rRNA loop -> translations -> gene sync -> sort.
  All earlier tests stay green (incl. KY753833 rRNA coordinates and the table2asn before/after test).

## Finding 14 restore the record on failure
- Test `test_failed_curation_leaves_record_unchanged` (main() in-process, Curator.run monkeypatched to mutate and raise; FAIL before: gene renamed, last feature dropped in the output).
- Fix: main() keeps `copy.deepcopy(rec)` before Curator(...).run() and puts it back into `recs` on exception (`import copy` moved to module level).

## Finding 20 sort key
- Test `test_origin_spanning_feature_sorts_by_its_larger_coordinate_part` (FAIL with the old `int(location.start)` key: wrapped trnY first).
- Fix: `Curator.sort_start`: for a compound feature with a part starting at 0, the sort position is the start of the part that
  ends at L (if there is one), else min part start; simple features unchanged. source/gene tie-breaks unchanged.

## Totals
- `python3 -m pytest -q tests/test_mito_curate.py tests/test_ncbi_validate.py`: 31 passed (tests/test_mito_curate.py 27 = 16 old + 11 new; tests/test_ncbi_validate.py 4 new).
  test_mito_rrna_hmm / test_mitom_cli / test_robustness: 53 passed (they import mito_curate).
- `ruff check mito_curate.py ncbi_validate.py tests/`: clean. mypy --ignore-missing-imports: 0 errors in mito_curate.py / ncbi_validate.py
  before and after (baseline HEAD copy and edited files both report only the pre-existing mito_env.py:255 error, not edited here).
- Every new test was run red first (finding 1 tbl + table2asn, 6, 7, 9, 10, 11, 12, 14, 15, 20); `test_rrna_trimmed_for_inside_gene_near_the_end`
  and the plus-strand tbl test are regression guards (green before and after).

## Effect on the 10 fixtures (tests/data, table2asn before -> after), compared with the HEAD scripts
- Same ERROR counts (9 -> 2, 7 types fixed, 0 new). WARNING 43 (HEAD) -> 45: NC_010217.1 orf345 (an intron ORF at 11945..12949 lying
  in cox1 intron 11548..13062; cox1 exon 2 11936..11944 is in the same frame) used to move its start 24 nt upstream (P1, overlapping that
  exon); with per-part spans the join always blocks, so the start is not moved and P5 marks it 5'-partial (TTA is no start codon):
  2 new PartialProblem warnings, no new error.
- NC_006902.1: R3 "rRNA overlaps a CDS by 1200 nt" is gone: that CDS is a join and the rRNA lies in its intron (hull overlap was spurious).

## Not done / notes
- Finding 6: `pcg_overlaps` only compares `simple()` CDS, so there is no hull there; unchanged. Partial (P5/P6) CDS also drop out of
  `simple()` and therefore out of P2 - pre-existing, not changed.
- Finding 11: start-codon / internal-stop checks on the joined CDS (reviewer's second suggestion) not added (not in the decision list).
- Finding 9: /transl_except other than TERM (none produced by this pipeline) is not honoured when translating.
