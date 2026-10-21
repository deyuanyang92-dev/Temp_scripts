# D fixes: mito_rescue.py / mito_intron.py notes (review C findings 2, 3, 4, 5, 15)

Files touched: mito_rescue.py, tests/test_mito_rescue.py, mito_intron.py (two note strings only),
tests/test_mito_intron.py. TDD: tests written first, each seen failing for the intended reason, then fixed.
No git commit.

## Finding 2 - containment release too weak (mito_rescue.rescue)
- Failing first: ND5 0..1800 + MITOS2 nad4 100..700 / nad2 800..1500 -> ['ND4','ND2'] (expected []);
  same on the opposite strand (nad6 / nad3); "only absent genes count" case returned ['CYTB'].
- Fix: new `MAX_LEN` table (nt, keyed by MITOS2 name; MitoZ names resolved through MITOZ_NAMES by `max_len()`),
  `CONTAINER_FACTOR = 1.5`. A CDS/rRNA releases its inside genes only if len > 1.5 x MAX_LEN of its own gene
  AND >= 2 distinct inside genes of other names that are absent from `have` lie wholly in it. Unknown own-gene
  name -> no release (conservative). tRNA never releases.
- Tests: test_short_container_never_released, test_container_release_counts_only_absent_genes,
  test_max_len_table_covers_all_genes; test_genes_inside_intron_spanning_cds_added adapted (container
  ND5 100..4900, SEQ now 5000 nt, second nad5 piece 4600..4900; other tests' coordinates unchanged and valid).
- Docstring rules paragraph rewritten.

## Finding 3 - `taken` not updated; tRNAs ignored
- Failing first: nad4 1000..1900 + nad5 1100..2200 -> ['ND4','ND5'] (expected ['ND5']); tRNA 1000..1070 +
  atp8 1000..1070 -> ['ATP8'] (expected []).
- Fix: tRNA positions added to `taken`; candidates processed longest-first, each accepted gene's positions
  appended to `taken` (50 % rule unchanged).
- Tests: test_overlapping_mitos_genes_longest_wins, test_trna_blocks_overlapping_prediction.

## Finding 5 - fragments joined regardless of distance
- Failing first: rrnL-a / rrnL-b became one rRNA feature; rrnL-a 1000..1200 + rrnL-b 1300..1700 gave
  (1000,1700) instead of the longest piece (1300,1700); CDS -a/-b 3200 nt apart or on opposite strands joined.
- Fix: rRNA lines are separate features (never joined). Non-rRNA fragments (CDS, also tRNA, same helper) are
  chained in letter order only when same strand and 0 <= gap <= MAX_FRAG_GAP (3000); otherwise separate
  features. Overlapping fragments (gap < 0) are also kept separate (a join with overlapping parts is invalid;
  this is my reading, beyond the literal "gap <= 3000").
- Tests: test_rrna_fragments_not_joined (3000..3160 instead of 5000..5160, same logic), 
  test_rrna_fragment_longest_piece_added, test_cds_fragments_joined_only_when_close_and_same_strand.
  A wrapped (origin-spanning) + strand gene was checked by script: location and translation correct.

## Finding 4 - /translation
- Failing first: BiopythonWarning (partial codon) with a 68-nt gene; internal-stop CDS still got a translation.
- Fix: `_translation()`: trim to whole codons, translate with the table, remove ONE trailing '*', first aa -> M
  when the first codon is in the table's start codons (Bio.Data.CodonTable), omit /translation when an internal
  '*' remains (a second trailing '*' counts as internal, so it is not hidden by rstrip).
- Tests: test_translation_trimmed_alt_start_is_m (table 5, ATT start, warnings as errors),
  test_translation_internal_stop_omitted, test_translation_non_start_first_codon_kept.
  test_missing_cds_added_mitoz_style adapted: its atp8 region was random sequence (internal stops), now a clean
  ORF, otherwise the new rule rightly omits the translation.

## Finding 15 - tool chatter in feature notes
- Failing first: note text and missing "hint" column.
- mito_rescue: NOTE = "annotated by similarity to a MITOS2 prediction (not found by the primary annotation)".
  The advice moved to a `hint` column of the report rows (HINT constant); rescue_report.tsv written by
  `main()` now has columns file, record, gene, location, status, hint. Failure rows have an empty hint.
- mito_intron: CDS note -> "exon-intron structure predicted by protein alignment (Exonerate, MFannot method)";
  intron note tail "mito_intron.py (MFannot method)" -> "MFannot method". No other change.
- Tests: test_note_is_neutral_and_hint_goes_to_report, test_report_has_hint_column,
  test_feature_notes_read_like_genbank_notes. The three selectors in test_mito_intron.py that found "our" CDS by
  `"mito_intron" in note` now compare with the new CDS note.

## Verification (final run)
- tests/test_mito_rescue.py: 21 passed. tests/test_mito_intron.py: 11 passed in 56 s (baseline before change: 10 passed,
  112 s). ruff check mito_rescue.py mito_intron.py tests/test_mito_rescue.py: clean.
  mypy --ignore-missing-imports mito_rescue.py: 0 errors in mito_rescue.py (mypy also follows the import into
  mito_env.py and reports 6 errors there, not my file).

## Not done / needs the owner of the other files
- batch_mitoz.py:3025-3026 writes the pipeline's rescue_report.tsv with an explicit 5-column header and key list;
  the new `hint` column does not appear there until "hint" is added to header and keys (extra row key is ignored,
  nothing breaks).
- benchmark/eval_cds.py:49 identifies mito_intron CDS by `"mito_intron" in note`; the new note has no script name.
  Change it to `"MFannot method" in n`.
- mito_intron.py:1015 label " - UNCONFIRMED: long intron from chained alignments, check manually" is still in
  the intron note (rare long-intron case; also the report's note column). Not changed: outside the stated
  decision, and it marks real uncertainty. ncbi_validate.INTERNAL_NOTE (another agent) matches "check .* manually"
  and would strip it from the .tbl.
