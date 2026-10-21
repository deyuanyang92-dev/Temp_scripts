"""Tests for mito_rrna_hmm.py (last-resort rRNA search with per-phylum profile HMMs).

Run: python3 -m pytest -q tests/test_mito_rrna_hmm.py   (the build/search test needs mafft, hmmbuild, nhmmer)
"""
import random
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_env  # noqa: E402
import mito_rrna_hmm as rh  # noqa: E402

sys.path.insert(0, str(HERE.parent / "scripts"))
import fetch_rrna_refs as fr  # noqa: E402

rnd = random.Random(5)


def rand(n):
    return "".join(rnd.choice("ACGT") for _ in range(n))


def rec(seq, feats=(), taxonomy=("Eukaryota", "Metazoa", "Ctenophora", "Tentaculata")):
    r = SeqRecord(Seq(seq), id="S1", description="")
    r.annotations.update(molecule_type="DNA", topology="circular", taxonomy=list(taxonomy))
    for kind, s, e, strand, gene in feats:
        loc = FeatureLocation(s, e, strand)
        r.features += [SeqFeature(loc, type="gene", qualifiers={"gene": [gene]}),
                       SeqFeature(loc, type=kind, qualifiers={"gene": [gene]})]
    return r


def hit(env_from, env_to, hmm_from, hmm_to, strand=1, score=100.0):
    return rh.Hit(env_from, env_to, strand, hmm_from, hmm_to, score)


def test_missing_rrnas_mitoz_and_ncbi_names():
    assert rh.missing(rec(rand(100), [("rRNA", 0, 50, 1, "s-rRNA")])) == ["rrnL"]
    assert rh.missing(rec(rand(100), [("rRNA", 0, 50, 1, "rrnL")])) == ["rrnS"]
    assert rh.missing(rec(rand(100), [("rRNA", 0, 40, 1, "rrnS"), ("rRNA", 50, 90, 1, "l-rRNA")])) == []


def test_gene_of_knows_the_aliases_of_mito_curate():
    """rrn12 / rrn16 / srRNA / lrRNA were unknown here although mito_curate (R3) recognises them: an existing
    feature named so was treated as missing and a second 16S was added."""
    def g(gene=None, product=None):
        q = {k: [v] for k, v in (("gene", gene), ("product", product)) if v}
        return rh.gene_of(SeqFeature(FeatureLocation(0, 10, 1), type="rRNA", qualifiers=q))
    assert [g(gene=x) for x in ("rrn12", "srRNA", "12S rRNA", "rrnS", "s-rRNA", "RNS")] == ["rrnS"] * 6
    assert [g(gene=x) for x in ("rrn16", "lrRNA", "16S rRNA", "rrnL", "l-rRNA", "RNL")] == ["rrnL"] * 6
    assert g(product="12S ribosomal RNA") == "rrnS" and g(product="16S ribosomal RNA") == "rrnL"
    assert g(product="small subunit ribosomal RNA") == "rrnS" and g(gene="foo") is None and g() is None
    # free text: "s-rrna" inside "16S-rRNA" is not the small subunit; "rnl" inside "trnL" is not an rRNA name
    assert g(product="mitochondrial 16S-rRNA") == "rrnL" and g(gene="16S-rRNA", product="large subunit ribosomal RNA") == "rrnL"
    assert g(gene="trnL") is None and g(gene="mt-rns") == "rrnS"
    assert rh.missing(rec(rand(100), [("rRNA", 0, 40, 1, "rrn12"), ("rRNA", 50, 90, 1, "lrRNA")])) == []


def test_phylum_from_taxonomy_or_override():
    lib = {"Ctenophora", "Nematoda"}
    assert rh.phylum_of(rec(rand(10)), lib) == "Ctenophora"
    assert rh.phylum_of(rec(rand(10), taxonomy=["Eukaryota", "Metazoa", "Porifera"]), lib) is None
    assert rh.phylum_of(rec(rand(10)), lib, override="Nematoda") == "Nematoda"


def test_fragments_of_one_gene_merged():
    """Two hits covering consecutive parts of the model (1-400, 420-900) in genome order are one gene."""
    best = rh.best_span([hit(1000, 1400, 1, 400), hit(1450, 1950, 420, 900, score=80)])
    assert (best.start, best.end) == (1000, 1950)


def test_duplicate_copies_not_merged():
    """User 2026-10-21: mitochondrial genes can be duplicated. Two hits over the same model columns are
    two copies, never one gene: only the better copy is kept."""
    best = rh.best_span([hit(1000, 1900, 1, 900, score=90), hit(2000, 2900, 1, 900, score=120)])
    assert (best.start, best.end) == (2000, 2900)


def test_model_order_must_match_genome_order():
    """Plus strand: the later genome piece must cover later model columns; minus strand the reverse."""
    assert (rh.best_span([hit(1000, 1400, 500, 900), hit(1450, 1850, 1, 400, score=80)]).end) == 1400
    m = rh.best_span([hit(1000, 1400, 500, 900, strand=-1), hit(1450, 1850, 1, 400, strand=-1, score=80)])
    assert (m.start, m.end) == (1000, 1850)


def test_far_pieces_not_merged():
    assert rh.best_span([hit(1000, 1400, 1, 400), hit(2000, 2500, 420, 900, score=80)]).end == 1400


def test_clusters_ranked_by_summed_score():
    """A gene in two weaker pieces (60 + 70 bits) beats one stray hit of 80 bits: clusters are compared by
    their summed score, not by their best single hit."""
    a1, a2, stray = hit(1000, 1300, 1, 300, score=60), hit(1350, 1930, 320, 900, score=70), hit(5000, 5200, 200, 400, score=80)
    best = rh.best_span([stray, a1, a2])
    assert (best.start, best.end, best.score) == (1000, 1930, 130)
    assert rh.best_span([a1, a2]) == best and rh.best_span([]) is None


def test_second_copy_is_the_runner_up_over_the_same_columns():
    c1, c2 = hit(1000, 1900, 1, 900, score=120), hit(2000, 2900, 1, 900, score=90)
    stray = hit(5000, 5200, 200, 400, score=80)                # a fragment, not a copy
    best, copy = rh.best_spans([c2, stray, c1])
    assert (best.start, copy.start, copy.end) == (1000, 2000, 2900)
    assert rh.best_spans([c1, stray])[1] is None and rh.best_spans([]) == (None, None)


def test_search_reports_second_copy_without_adding_it(tmp_path):
    lib = _lib(tmp_path / "lib")
    r = rec(rand(5000), [("CDS", 0, 800, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")])
    rows = rh.search_record(r, lib, run=lambda hmm, fa: [hit(1030, 1910, 1, 880, score=199), hit(2746, 3626, 1, 880, score=175)])
    assert [(x["gene"], x["status"]) for x in rows] == [("l-rRNA", "added"),
                                                         ("l-rRNA", "second copy at 2747..3626 not added")]
    assert sum(f.type == "rRNA" and f.qualifiers["gene"] == ["l-rRNA"] for f in r.features) == 1


def _lib(tmp_path, median=900, leng=880, p5=800, p95=1000):
    """Fake library: manifest row plus a stub .hmm that only has the header line LENG (read by the module).
    p5 / p95 None = a manifest from before these columns existed."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    cols, vals = ("\tp5_len\tp95_len", f"\t{p5}\t{p95}") if p5 is not None else ("", "")
    (tmp_path / "manifest.tsv").write_text(f"phylum\tgene\tn_ref\tn_used\tmedian_len{cols}\n"
                                           f"Ctenophora\trrnL\t30\t30\t{median}{vals}\n")
    (tmp_path / "Ctenophora.rrnL.hmm").write_text(f"HMMER3/f [3.4]\nNAME  Ctenophora.rrnL\nLENG  {leng}\nALPH  DNA\n")
    return tmp_path


def test_search_adds_rrna_with_note(tmp_path):
    lib = _lib(tmp_path / "lib")
    r = rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")])
    rows = rh.search_record(r, lib, run=lambda hmm, fa: [hit(2000, 2880, 1, 880)])
    assert [(x["gene"], x["status"]) for x in rows] == [("l-rRNA", "added")]
    f = [f for f in r.features if f.type == "rRNA" and f.qualifiers["gene"] == ["l-rRNA"]][0]
    assert (int(f.location.start), int(f.location.end)) == (2000, 2880)
    assert f.qualifiers["product"] == ["16S ribosomal RNA"] and "profile HMM" in f.qualifiers["note"][0]
    assert any(g.type == "gene" and g.qualifiers["gene"] == ["l-rRNA"] for g in r.features)


def test_feature_note_is_submission_style_and_the_to_do_stays_in_the_report(tmp_path):
    """The /note ends up in the NCBI .tbl and the public record: no tool name, no instruction to the submitter."""
    lib = _lib(tmp_path / "lib")
    r = rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")])
    rows = rh.search_record(r, lib, run=lambda hmm, fa: [hit(2000, 2880, 1, 880)])
    for f in r.features:
        if f.qualifiers.get("gene") == ["l-rRNA"]:
            assert f.qualifiers.get("note", [None])[0] in (None, "annotated by profile HMM similarity "
                                                                  "(per-phylum rRNA model, Ctenophora, 30 references)")
    (f,) = _added_rrna(r)
    assert "check" not in f.qualifiers["note"][0] and "mito_rrna_hmm" not in f.qualifiers["note"][0]
    assert "check boundaries manually" in rows[0]["note"]


def test_search_rejects_hit_on_existing_gene_short_or_long(tmp_path):
    lib = _lib(tmp_path / "lib")
    for h, why in ((hit(200, 1300, 1, 880), "existing"), (hit(2000, 2150, 1, 880), "short"),
                   (hit(1600, 3200, 1, 880), "long")):
        r = rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")])
        rows = rh.search_record(r, lib, run=lambda hmm, fa, h=h: [h])
        assert rows[0]["status"].startswith("rejected") and why in rows[0]["status"], (why, rows)
        assert not any(f.qualifiers.get("gene") == ["l-rRNA"] for f in r.features)


def _added_rrna(r):
    return [f for f in r.features if f.type == "rRNA" and f.qualifiers["gene"] == ["l-rRNA"]]


def test_low_model_coverage_is_rejected(tmp_path):
    """NC_016117.1 case of the review: a 336-nt hit over 38 % of the model was added as a complete 16S."""
    lib = _lib(tmp_path / "lib")
    r = rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")])
    rows = rh.search_record(r, lib, run=lambda hmm, fa: [hit(2000, 2336, 1, 336)])
    assert rows[0]["status"].startswith("rejected: low coverage") and rows[0]["coverage"] == "0.38", rows
    assert not _added_rrna(r)


@pytest.mark.parametrize("strand,cols,status,start_cls,end_cls", [
    (1, (261, 880), "added (partial 5')", BeforePosition, type(None)),
    (-1, (261, 880), "added (partial 5')", type(None), AfterPosition),
    (1, (1, 620), "added (partial 3')", type(None), AfterPosition),
    (-1, (1, 620), "added (partial 3')", BeforePosition, type(None)),
    (1, (201, 700), "added (partial 5' and 3')", BeforePosition, AfterPosition),
    (-1, (201, 700), "added (partial 5' and 3')", BeforePosition, AfterPosition)])
def test_partial_ends_are_marked_on_the_correct_side(tmp_path, strand, cols, status, start_cls, end_cls):
    """hmm_from > 30 = the 5' end of the gene is missing; hmm_to < LENG - 30 = the 3' end. 5' is the genome
    start on the plus strand and the genome end on the minus strand."""
    lib = _lib(tmp_path / "lib")
    r = rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")])
    rows = rh.search_record(r, lib, run=lambda hmm, fa: [hit(2000, 2000 + cols[1] - cols[0] + 1, cols[0], cols[1], strand)])
    assert rows[0]["status"] == status, rows
    for f in _added_rrna(r) + [g for g in r.features if g.type == "gene" and g.qualifiers["gene"] == ["l-rRNA"]]:
        for pos, cls in ((f.location.start, start_cls), (f.location.end, end_cls)):
            assert (type(pos) is cls) if cls is not type(None) else not isinstance(pos, (BeforePosition, AfterPosition))


def test_full_length_hit_has_no_partial_mark_and_reports_coverage(tmp_path):
    lib = _lib(tmp_path / "lib")
    r = rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")])
    rows = rh.search_record(r, lib, run=lambda hmm, fa: [hit(2000, 2870, 11, 870)])    # within 30 columns of both ends
    assert (rows[0]["status"], rows[0]["coverage"]) == ("added", "0.98")


def test_unreadable_model_is_a_failed_row(tmp_path):
    lib = _lib(tmp_path / "lib")
    (lib / "Ctenophora.rrnL.hmm").write_text("x")
    r = rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")])
    rows = rh.search_record(r, lib, run=lambda hmm, fa: [hit(2000, 2880, 1, 880)])
    assert rows[0]["status"].startswith("failed: ") and "LENG" in rows[0]["status"], rows


def _status(tmp_path, h, **lib_kw):
    r = rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")])
    return rh.search_record(r, _lib(tmp_path / "lib", **lib_kw), run=lambda hmm, fa: [h])[0]["status"]


def test_length_window_from_p5_and_p95_not_the_median(tmp_path):
    """Upper limit 1.3 x p95 (here 1300), lower limit 0.7 x p5 (560) for spans with both ends complete."""
    assert _status(tmp_path, hit(2000, 3290, 1, 880)) == "added"                       # 1290 nt
    assert _status(tmp_path, hit(2000, 3320, 1, 880)).startswith("rejected: long")     # 1320 nt
    assert _status(tmp_path, hit(2000, 2570, 1, 880)) == "added"                       # 570 nt
    assert _status(tmp_path, hit(2000, 2500, 1, 880)).startswith("rejected: short")    # 500 nt
    # a long gene is no outlier when the references say so: p95 1200 -> limit 1560 (the old 1.5 x median 900 = 1350)
    assert _status(tmp_path, hit(2000, 3500, 1, 880), p95=1200) == "added"


def test_partial_span_is_not_held_to_the_complete_gene_lower_limit(tmp_path):
    """A 400-nt span covering 55 % of the model with the 5' end missing is a partial gene, not 'short'."""
    assert _status(tmp_path, hit(2000, 2400, 401, 880)) == "added (partial 5')"
    assert _status(tmp_path, hit(2000, 2150, 401, 880)).startswith("rejected: short")  # but < 200 nt is


def test_manifest_without_p5_p95_keeps_the_old_median_rule(tmp_path):
    assert _status(tmp_path, hit(2000, 3340, 1, 880), p5=None, p95=None) == "added"      # 1340 <= 1.5 x 900
    assert _status(tmp_path, hit(2000, 3360, 1, 880), p5=None, p95=None).startswith("rejected: long")


def _refs_gb(path, lens, gene_type="rRNA"):
    recs = []
    for i, n in enumerate(lens):
        r = SeqRecord(Seq(rand(400 + n + 400)), id=f"R{i}.1", description="")
        r.annotations.update(molecule_type="DNA", topology="circular", taxonomy=["Eukaryota", "Metazoa", "Ctenophora"])
        r.features.append(SeqFeature(FeatureLocation(400, 400 + n, 1), type=gene_type,
                                     qualifiers={"product": ["16S ribosomal RNA"]}))
        recs.append(r)
    SeqIO.write(recs, str(path), "genbank")


def test_manifest_only_rebuild_keeps_models_and_adds_percentiles(tmp_path):
    """build(manifest_only=True) recomputes manifest.tsv (n_ref, median, p5, p95) from the references without
    mafft / hmmbuild and without touching the .hmm files; reference rRNAs up to 4000 nt are used."""
    lib = tmp_path / "lib"
    lib.mkdir()
    stub = "HMMER3/f [3.4]\nNAME  Ctenophora.rrnL\nLENG  880\nNSEQ  4\nALPH  DNA\n"
    (lib / "Ctenophora.rrnL.hmm").write_text(stub)
    (lib / "manifest.tsv").write_text("phylum\tgene\tn_ref\tn_used\tmedian_len\n"
                                      "Ctenophora\trrnL\t4\t4\t650\nNematoda\trrnL\t9\t9\t1000\n")
    gb = tmp_path / "Ctenophora.gb"
    _refs_gb(gb, [500, 600, 700, 800, 3500, 4500])                  # 4500 > cap 4000 is left out
    rows = rh.build([gb], lib, manifest_only=True)
    assert stub == (lib / "Ctenophora.rrnL.hmm").read_text() and not (lib / "Ctenophora.rrnS.hmm").exists()
    assert [(r["phylum"], r["gene"], r["n_ref"], r["n_used"], r["median_len"], r["p5_len"], r["p95_len"])
            for r in rows] == [("Ctenophora", "rrnL", 5, 4, 700, 520, 2960)]
    man = rh.read_manifest(lib)
    assert man[("Ctenophora", "rrnL")]["p95_len"] == "2960" and man[("Nematoda", "rrnL")]["n_ref"] == "9"


L = 5000                                  # genome length in the origin tests; the gene is 975 nt, 437 before the origin


def _wrap(tmp_path, hits, feats=(("CDS", 1800, 2500, 1, "COX1"), ("rRNA", 3000, 3600, 1, "s-rRNA")), topology="circular"):
    """Search a circular L-nt record with canned hits; returns (record, rows, the FASTA that was searched)."""
    lib = _lib(tmp_path / "lib", median=975, leng=975, p5=900, p95=1000)
    r = rec(rand(L), list(feats))
    r.annotations["topology"] = topology
    seen = []

    def run(hmm, fa):
        seen.append(fa.read_text())
        return hits
    return r, rh.search_record(r, lib, run=run), seen[0]


def test_circular_record_is_searched_with_its_first_2000_nt_appended(tmp_path):
    r, _, fa = _wrap(tmp_path, [])
    assert fa.split("\n")[1] == str(r.seq) + str(r.seq[:2000])
    r, _, fa = _wrap(tmp_path / "lin", [], topology="linear")
    assert fa.split("\n")[1] == str(r.seq)


def test_gene_across_the_origin_becomes_a_join_location(tmp_path):
    r, rows, _ = _wrap(tmp_path, [hit(L - 437, L + 538, 1, 975)])
    assert [x["status"] for x in rows] == ["added"] and rows[0]["location"] == "4564..5000,1..538+"
    (f,) = _added_rrna(r)
    assert isinstance(f.location, CompoundLocation)
    assert [(int(p.start), int(p.end), p.strand) for p in f.location.parts] == [(L - 437, L, 1), (0, 538, 1)]


def test_origin_gene_on_the_minus_strand_lists_parts_in_transcript_order(tmp_path):
    r, rows, _ = _wrap(tmp_path, [hit(L - 437, L + 538, 1, 975, strand=-1)])
    assert rows[0]["status"] == "added"
    (f,) = _added_rrna(r)
    assert [(int(p.start), int(p.end), p.strand) for p in f.location.parts] == [(0, 538, -1), (L - 437, L, -1)]
    assert len(f.extract(r.seq)) == 975


def test_pieces_across_the_origin_join_and_the_duplicate_in_the_first_part_is_not_a_copy(tmp_path):
    """nhmmer on genome+2000 nt can give the gene as two pieces, the second starting at the origin copy
    (>= L). Together with the real hit of that second piece (0..538) there must be ONE gene, no 'second copy'."""
    hits = [hit(L - 437, L, 1, 437), hit(L, L + 538, 438, 975), hit(0, 538, 438, 975)]
    r, rows, _ = _wrap(tmp_path, hits)
    assert [x["status"] for x in rows] == ["added"]
    assert [(int(p.start), int(p.end)) for p in _added_rrna(r)[0].location.parts] == [(L - 437, L), (0, 538)]


def test_hits_only_in_the_appended_copy_are_duplicates_and_dropped(tmp_path):
    """A gene inside the first 2000 nt is found twice (real + appended copy): one row, not a second copy."""
    r, rows, _ = _wrap(tmp_path, [hit(500, 1475, 1, 975, score=100), hit(L + 500, L + 1475, 1, 975, score=100)])
    assert [(x["status"], x["location"]) for x in rows] == [("added", "501..1475+")]
    assert rh.best_spans([hit(L + 500, L + 1475, 1, 975)], L) == (None, None)


def test_overlap_with_existing_genes_is_tested_modulo_the_genome_length(tmp_path):
    r, rows, _ = _wrap(tmp_path, [hit(L - 437, L + 538, 1, 975)], feats=(("CDS", 0, 520, 1, "COX1"), ("rRNA", 3000, 3600, 1, "s-rRNA")))
    assert rows[0]["status"] == "rejected: lies on an existing gene" and not _added_rrna(r)


def test_partial_origin_gene_marks_the_outer_ends(tmp_path):
    r, rows, _ = _wrap(tmp_path, [hit(L - 437, L + 300, 261, 975)])                       # 5' end missing, + strand
    assert rows[0]["status"] == "added (partial 5')"
    a, b = _added_rrna(r)[0].location.parts
    assert isinstance(a.start, BeforePosition) and not isinstance(b.end, AfterPosition)
    r, rows, _ = _wrap(tmp_path / "m", [hit(L - 437, L + 300, 261, 975, strand=-1)])      # minus: 5' end is the genome end
    a, b = _added_rrna(r)[0].location.parts
    assert isinstance(a.end, AfterPosition) and int(a.end) == 300 and not isinstance(b.start, BeforePosition)


def test_report_file_has_coverage_and_note_columns(tmp_path):
    lib = _lib(tmp_path / "lib")
    d = tmp_path / "final"
    _gbf(d / "S1.gbf", rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")]))
    rows = rh.search_dir(d, lib, run=lambda hmm, fa: [hit(2000, 2620, 261, 880)])
    rh.write_report(rows, d / "rrna_hmm_report.tsv")
    head, line = (d / "rrna_hmm_report.tsv").read_text().splitlines()[:2]
    assert head.split("\t") == ["file", "record", "phylum", "gene", "location", "status", "coverage", "note"]
    assert line.split("\t")[5:7] == ["added (partial 5')", "0.70"] and "check boundaries manually" in line


def test_search_skips_when_present_or_no_library(tmp_path):
    lib = _lib(tmp_path / "lib")
    never = lambda hmm, fa: (_ for _ in ()).throw(AssertionError("nhmmer must not run"))  # noqa: E731
    full = rec(rand(3000), [("rRNA", 0, 800, 1, "s-rRNA"), ("rRNA", 900, 1900, 1, "l-rRNA")])
    assert rh.search_record(full, lib, run=never) == []
    other = rec(rand(3000), [("rRNA", 0, 800, 1, "s-rRNA")], taxonomy=["Eukaryota", "Metazoa", "Porifera"])
    rows = rh.search_record(other, lib, run=never)
    assert rows[0]["status"] == "skipped: no library for this phylum"


def test_search_dir_in_place(tmp_path):
    lib = _lib(tmp_path / "lib")
    d = tmp_path / "final"
    d.mkdir()
    SeqIO.write(rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")]),
                str(d / "S1.gbf"), "genbank")
    (d / "S1.gbf").chmod(0o644)
    rows = rh.search_dir(d, lib, run=lambda hmm, fa: [hit(2000, 2880, 1, 880)])
    assert [(x["record"], x["gene"], x["status"]) for x in rows] == [("S1", "l-rRNA", "added")]
    genes = [f.qualifiers["gene"][0] for f in SeqIO.read(str(d / "S1.gbf"), "genbank").features if f.type == "rRNA"]
    assert sorted(genes) == ["l-rRNA", "s-rRNA"]
    assert (d / "S1.gbf").stat().st_mode & 0o777 == 0o644          # atomic rewrite keeps the file mode


def _fake_nhmmer(tmp_path, monkeypatch, body):
    exe = tmp_path / "nhmmer"
    exe.write_text("#!/bin/sh\n" + body + "\n")
    exe.chmod(0o755)
    monkeypatch.setattr(rh.mito_env, "find", lambda name: str(exe))


def test_nhmmer_error_includes_stderr(tmp_path, monkeypatch):
    """A failing nhmmer must say why (its stderr), not only 'non-zero exit status 1'."""
    _fake_nhmmer(tmp_path, monkeypatch, "echo 'Error reading query file (7)' >&2; exit 1")
    fa = tmp_path / "q.fa"
    fa.write_text(">x\nACGT\n")
    with pytest.raises(RuntimeError, match=r"Error reading query file \(7\)"):
        rh.nhmmer(tmp_path / "m.hmm", fa)


def test_search_record_failure_becomes_report_row(tmp_path):
    lib = _lib(tmp_path / "lib")

    def boom(hmm, fa):
        raise RuntimeError("nhmmer exploded: bad model")

    def boom_cpe(hmm, fa):
        raise subprocess.CalledProcessError(1, "nhmmer", stderr="Error reading query file")

    for run, text in ((boom, "nhmmer exploded: bad model"), (boom_cpe, "Error reading query file")):
        r = rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")])
        rows = rh.search_record(r, lib, run=run)
        assert len(rows) == 1 and rows[0]["gene"] == "l-rRNA"
        assert rows[0]["status"].startswith("failed: ") and text in rows[0]["status"], rows


def _gbf(path, r):
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(r, str(path), "genbank")


def test_search_dir_survives_a_failing_record_and_reads_subdirectories(tmp_path):
    """One bad record must not stop the others; records in 04.final_gb/no_cox1_gene/ are searched too."""
    lib = _lib(tmp_path / "lib")
    d = tmp_path / "final"
    ok = rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")])
    bad = rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")])
    bad.id = "BAD"
    _gbf(d / "a.gbf", bad)
    _gbf(d / "no_cox1_gene" / "S1.gbf", ok)
    _gbf(d / "no_cox1_gene" / "S2.gb", dict_id(ok, "S2"))

    def run(hmm, fa):
        if fa.read_text().startswith(">BAD"):
            raise RuntimeError("boom")
        return [hit(2000, 2880, 1, 880)]

    rows = rh.search_dir(d, lib, run=run)
    got = {(x["file"], x["record"]): x["status"] for x in rows}
    assert got[("a.gbf", "BAD")] == "failed: boom"
    assert got[("no_cox1_gene/S1.gbf", "S1")] == "added" and got[("no_cox1_gene/S2.gb", "S2")] == "added"
    assert {f.qualifiers["gene"][0] for f in SeqIO.read(str(d / "no_cox1_gene" / "S1.gbf"), "genbank").features
            if f.type == "rRNA"} == {"l-rRNA", "s-rRNA"}


def dict_id(r, new_id):
    r2 = r[:]
    r2.id, r2.name = new_id, new_id
    r2.annotations = dict(r.annotations)
    return r2


def test_search_dir_failure_outside_search_record_is_a_row(tmp_path, monkeypatch):
    lib = _lib(tmp_path / "lib")
    d = tmp_path / "final"
    _gbf(d / "a.gbf", rec(rand(5000), [("CDS", 0, 1500, 1, "COX1")]))

    def boom(*a, **k):
        raise RuntimeError("library unreadable")

    monkeypatch.setattr(rh, "search_record", boom)
    rows = rh.search_dir(d, lib)
    assert [(x["file"], x["record"], x["status"]) for x in rows] == [("a.gbf", "S1", "failed: library unreadable")]


def test_gbf_written_atomically(tmp_path, monkeypatch):
    """A failure while writing the GenBank file must leave the original untouched (and no temp file behind)."""
    lib = _lib(tmp_path / "lib")
    d = tmp_path / "final"
    _gbf(d / "a.gbf", rec(rand(5000), [("CDS", 0, 1500, 1, "COX1"), ("rRNA", 4000, 4600, 1, "s-rRNA")]))
    before = (d / "a.gbf").read_bytes()
    real = SeqIO.write

    def broken(recs, dest, fmt):
        if hasattr(dest, "write"):
            dest.write("LOCUS partial")
        else:
            Path(dest).write_text("LOCUS partial")
        raise ValueError("cannot serialise")

    monkeypatch.setattr(rh.SeqIO, "write", broken)
    rows = rh.search_dir(d, lib, run=lambda hmm, fa: [hit(2000, 2880, 1, 880)])
    monkeypatch.setattr(rh.SeqIO, "write", real)
    assert (d / "a.gbf").read_bytes() == before
    assert sorted(p.name for p in d.iterdir()) == ["a.gbf"]
    assert rows and rows[0]["status"].startswith("failed: ")


def _mutate(s, rate):
    return "".join(rnd.choice("ACGT".replace(c, "")) if rnd.random() < rate else c for c in s)


@pytest.mark.skipif(not all(mito_env.find(t) for t in ("mafft", "hmmbuild", "nhmmer")),
                    reason="mafft / hmmbuild / nhmmer not found")
def test_build_and_search_real_tools(tmp_path):
    """End to end: library from 6 reference genomes of one phylum, then a new genome whose 16S differs by
    ~15 % from every reference is found by nhmmer."""
    anc = rand(900)
    refs = []
    for i in range(6):
        r = SeqRecord(Seq(rand(400) + _mutate(anc, 0.12) + rand(400)), id=f"R{i}", description="")
        r.annotations.update(molecule_type="DNA", topology="circular", organism=f"Genus{i} sp{i}",
                             taxonomy=["Eukaryota", "Metazoa", "Ctenophora"])
        r.features.append(SeqFeature(FeatureLocation(400, 400 + len(r.seq) - 800, 1), type="rRNA",
                                     qualifiers={"product": ["16S ribosomal RNA"]}))
        refs.append(r)
    gb = tmp_path / "Ctenophora.gb"
    SeqIO.write(refs, str(gb), "genbank")
    lib = tmp_path / "lib"
    rh.build([gb], lib)
    assert (lib / "Ctenophora.rrnL.hmm").exists()
    q = rec(rand(1500) + _mutate(anc, 0.15) + rand(1500), [("CDS", 0, 1400, 1, "COX1")])
    rows = rh.search_record(q, lib)
    added = [x for x in rows if x["gene"] == "l-rRNA"]
    assert added and added[0]["status"] == "added", rows
    f = [f for f in q.features if f.type == "rRNA"][0]
    assert abs(int(f.location.start) - 1500) < 60 and abs(int(f.location.end) - 2400) < 60


# ---- scripts/fetch_rrna_refs.py (reference download; the NCBI calls are replaced by a fake `get`) ----

def _gb_text(acc, genus="Genus"):
    r = SeqRecord(Seq(rand(60)), id=acc, name=acc.split(".")[0], description="mitochondrion, complete genome")
    r.annotations.update(molecule_type="DNA", topology="circular", organism=f"{genus} sp", taxonomy=["Eukaryota", "Metazoa"])
    import io
    h = io.StringIO()
    SeqIO.write(r, h, "genbank")
    return h.getvalue()


def _fake_ncbi(records):
    """records: {uid: (accession.version, genus)} -> a get(url) like the E-utilities."""
    import json

    def get(url):
        if "esearch" in url:
            return json.dumps({"esearchresult": {"count": str(len(records)), "idlist": list(records)}})
        ids = url.split("&id=")[1].split("&")[0].split(",")
        return "".join(_gb_text(*records[i]) for i in ids)
    return get


def test_fetch_exclusions_from_the_benchmark_set(tmp_path):
    tsv = tmp_path / "set.tsv"
    tsv.write_text("accession\tsource_accession\tphylum\tgenus\nNC_016754.1\tJN710452\tAcanthocephala\tOncicola\n"
                   "NC_022921.1\t\tAcanthocephala\tPallisentis\n")
    acc, gen = fr.load_set(tsv)
    assert acc == {"NC_016754", "JN710452", "NC_022921"} and gen == {"Oncicola", "Pallisentis"}


def test_fetch_phylum_excludes_samples_and_is_reproducible():
    recs = {str(i): (f"AB{i:06d}.1", "Oncicola" if i == 3 else "Genus") for i in range(1, 21)}
    get = _fake_ncbi(recs)
    texts, info = fr.fetch_phylum("Foo", get, per=5, seed=1, exclude_acc={"AB000001"}, pause=0)
    accs = [fr.accession_of(x) for x in texts]
    assert len(accs) == 5 and "AB000001.1" not in accs and info["n_hits"] == 20 and info["n_kept"] == 5
    assert accs == [fr.accession_of(x) for x in fr.fetch_phylum("Foo", get, per=5, seed=1, exclude_acc={"AB000001"}, pause=0)[0]]
    allgen, _ = fr.fetch_phylum("Foo", get, per=50, seed=1, exclude_gen={"Oncicola"}, pause=0)
    assert len(allgen) == 19 and "AB000003.1" not in [fr.accession_of(x) for x in allgen]
    assert "mitochondrion[Filter]" in info["query"] and info["query"].startswith("Foo[Organism]")


def test_fetch_main_writes_gb_accession_list_and_log(tmp_path, monkeypatch):
    recs = {str(i): (f"AB{i:06d}.1", "Genus") for i in range(1, 8)}
    monkeypatch.setattr(fr, "ncbi_get", _fake_ncbi(recs))
    out = tmp_path / "refs"
    assert fr.main(["-o", str(out), "--phyla", "Foo", "Bar", "--per-phylum", "4", "--seed", "3", "--pause", "0"]) == 0
    for p in ("Foo", "Bar"):
        accs = (out / f"refs_{p}.txt").read_text().split()
        assert len(accs) == 4 and accs == [r.id for r in SeqIO.parse(str(out / f"{p}.gb"), "genbank")]
    log = (out / "fetch_log.tsv").read_text().splitlines()
    assert log[0].split("\t")[:3] == ["phylum", "query", "date"] and len(log) == 3


def test_fetch_metazoan_phyla_list_matches_batch_mitoz():
    """The script keeps its own copy (it must not import batch_mitoz); the copy must not go stale."""
    import batch_mitoz
    assert set(fr.METAZOAN_PHYLA) == set(batch_mitoz.METAZOAN_PHYLA)
