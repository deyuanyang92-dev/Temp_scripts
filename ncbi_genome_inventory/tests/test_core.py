"""Offline tests. Fixtures are trimmed live NCBI responses for Acholoe squamosa (Polynoidae, 2026-10-28)."""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import core  # noqa: E402
import genome_inventory as gi  # noqa: E402

FIX = HERE / "fixtures"


def load(name):
    return json.loads((FIX / name).read_text(encoding="utf-8"))


def run(**kw):
    base = {"Run": "R1", "LibraryStrategy": "WGS", "LibrarySource": "GENOMIC", "LibrarySelection": "RANDOM",
            "Platform": "ILLUMINA", "Model": "Illumina NovaSeq 6000", "LibraryName": ""}
    base.update(kw)
    return base


# ---------------------------------------------------------------- classification

def test_hic_by_strategy_and_by_keyword():
    assert core.classify_run(run(LibraryStrategy="Hi-C"))[0] == "Hi-C"
    assert core.classify_run(run(LibraryStrategy="OTHER"), {"library_protocol": "Hi-C - Arima v2"})[0] == "Hi-C"
    assert core.classify_run(run(), {"experiment_title": "Omni-C library"})[0] == "Hi-C"


def test_hic_keyword_ignored_for_rna():
    # an RNA run whose title mentions Arima must stay RNA-seq
    r = run(LibraryStrategy="RNA-Seq", LibrarySource="TRANSCRIPTOMIC")
    assert core.classify_run(r, {"experiment_title": "arima"})[0] == "RNA-seq"


def test_pacbio_rules():
    pb = dict(Platform="PACBIO_SMRT")
    assert core.classify_run(run(Model="Sequel II", **pb), {"library_protocol": "PacBio - HiFi"})[0] == "HiFi"
    assert core.classify_run(run(Model="Sequel", **pb), {"experiment_title": "subreads"})[0] == "CLR"
    assert core.classify_run(run(Model="Revio", **pb))[0] == "HiFi"
    assert core.classify_run(run(Model="PacBio RS II", **pb))[0] == "CLR"
    assert core.classify_run(run(Model="Sequel IIe", **pb))[0] == "HiFi"          # inferred, rule says so
    assert "inferred" in core.classify_run(run(Model="Sequel IIe", **pb))[1]
    assert core.classify_run(run(Model="Sequel II", **pb))[0] == "PacBio_unknown"
    assert core.classify_run(run(Model="Sequel II", LibrarySource="TRANSCRIPTOMIC", **pb))[0] == "Iso-Seq"


def test_ont_rna_shortwgs_10x_metagenomic():
    assert core.classify_run(run(Platform="OXFORD_NANOPORE"))[0] == "ONT"
    assert core.classify_run(run(Platform="OXFORD_NANOPORE", LibrarySource="TRANSCRIPTOMIC"))[0] == "ONT_RNA"
    assert core.classify_run(run(LibraryStrategy="RNA-Seq", LibrarySource="TRANSCRIPTOMIC"))[0] == "RNA-seq"
    assert core.classify_run(run())[0] == "short_read_WGS"
    assert core.classify_run(run(Platform="BGISEQ"))[0] == "short_read_WGS"
    assert core.classify_run(run(), {"library_protocol": "10X Genomics Chromium"})[0] == "10x_linked"
    assert core.classify_run(run(LibrarySource="METAGENOMIC"))[0] == "Other"
    assert core.classify_run(run(LibraryStrategy="AMPLICON"))[0] == "Other"


# ---------------------------------------------------------------- NCBI columns / assemblies

def test_ncbi_columns_match_user_export_names():
    reps = load("acholoe_assemblies.json")
    row = core.ncbi_columns(reps[0])
    assert list(row) == core.NCBI_COLS and len(core.NCBI_COLS) == 41
    prim = core.ncbi_columns(next(r for r in reps if r["accession"] == "GCA_949317995.1"))
    assert prim["Assembly Stats Total Sequence Length"] == 1481470071          # numeric, not text
    assert prim["Assembly BioProject Accession"] == "PRJEB60118"
    assert prim["Assembly Sequencing Tech"] == "PacBio,Arima2"


def test_merge_copies_gcf_annotation_onto_gca():
    gca = {"accession": "GCA_1.1", "paired_accession": "GCF_1.1", "organism": {"tax_id": 1},
           "assembly_info": {}, "assembly_stats": {}}
    gcf = {"accession": "GCF_1.1", "paired_accession": "GCA_1.1", "organism": {"tax_id": 1},
           "assembly_info": {}, "assembly_stats": {},
           "annotation_info": {"name": "NCBI Annotation Release 100",
                               "stats": {"gene_counts": {"total": 20000}}}}
    rows = core.merge_gca_gcf([gca, gcf])
    assert len(rows) == 1
    assert rows[0]["Annotation Name"] == "NCBI Annotation Release 100"
    assert rows[0]["Annotation Count Gene Total"] == 20000
    assert rows[0]["Assembly Paired Assembly Accession"] == "GCF_1.1"


def test_representative_prefers_primary_over_alternate():
    rows = core.merge_gca_gcf(load("acholoe_assemblies.json"))
    assert core.representative(rows)["assembly_accession"] == "GCA_949317995.1"


# ---------------------------------------------------------------- linking (the Sanger/DToL case)

def test_sanger_raw_data_linked_through_sibling_project():
    """HiFi/Hi-C runs sit in PRJEB59385 with other BioSamples; the assembly is PRJEB60118. Both hang under the
    species umbrella PRJEB59386 -> must link as sibling_project (BioSample alone would miss them)."""
    asm = next(r for r in core.merge_gca_gcf(load("acholoe_assemblies.json"))
               if r["assembly_accession"] == "GCA_949317995.1")
    sra, bp = load("acholoe_sra.json"), load("acholoe_bioprojects.json")
    ev = {r["Run"]: core.link_evidence(r, asm, bp) for r in sra["runinfo"]}
    assert ev == {"ERR10841326": "sibling_project", "ERR10841327": "sibling_project",
                  "ERR10851523": "sibling_project"}
    cats = {r["Run"]: core.classify_run(r, sra["xml"][r["Run"]])[0] for r in sra["runinfo"]}
    assert cats == {"ERR10841326": "HiFi", "ERR10841327": "HiFi", "ERR10851523": "Hi-C"}


def test_unrelated_project_not_linked():
    asm = {"biosample": "SAMX", "bioproject": "PRJ_A", "bioproject_parents": "PRJ_P"}
    assert core.link_evidence({"BioSample": "SAMY", "BioProject": "PRJ_B"}, asm, {"PRJ_B": ["PRJ_Q"]}) == ""
    assert core.link_evidence({"BioSample": "SAMX", "BioProject": "PRJ_B"}, asm, {}) == "biosample"
    assert core.link_evidence({"BioSample": "", "BioProject": "PRJ_P"}, asm, {}) == "bioproject"


def test_parse_bioproject_parents():
    xml = """<RecordSet><DocumentSummary><Project><ProjectID><ArchiveID accession="PRJEB59385" archive="EBI"
    id="929871"/></ProjectID></Project><ProjectLinks><Link><ProjectIDRef archive="EBI" id="929871"
    accession="PRJEB59385"/><Hierarchical type="TopAdmin"><MemberID archive="EBI" id="929872"
    accession="PRJEB59386"/></Hierarchical></Link></ProjectLinks></DocumentSummary></RecordSet>"""
    assert core.parse_bioproject_parents(xml) == {"PRJEB59385": ["PRJEB59386"]}


def test_mito_from_sequence_report():
    reps = [{"genbank_accession": "OX439053.1", "length": 15641, "assigned_molecule_location_type": "Mitochondrion"},
            {"genbank_accession": "OX439035.1", "length": 1, "assigned_molecule_location_type": "Chromosome"}]
    assert core.mito_from_seqreport(reps) == [("OX439053.1", 15641)]


def test_species_of_handles_subspecies():
    tax = {"222": {"ancestors": ["1", "111"]}}
    assert gi.species_of("222", ["111", "333"], tax) == ["111"]
    assert gi.species_of("333", ["111", "333"], tax) == ["333"]


def test_mitogenome_only_assembly_flagged_and_not_representative():
    def rep(acc, length, level):
        return {"accession": acc, "organism": {"tax_id": 1}, "assembly_info": {"assembly_level": level},
                "assembly_stats": {"total_sequence_length": str(length), "contig_n50": length}}
    rows = core.merge_gca_gcf([rep("GCA_965119645.1", 18801, "Complete Genome"),
                               rep("GCA_903813345.2", 500_000_000, "Scaffold")])
    flags = {r["assembly_accession"]: r["assembly_content_derived"] for r in rows}
    assert flags == {"GCA_965119645.1": "mitogenome-only (< 1 Mb)", "GCA_903813345.2": "nuclear genome"}
    assert core.representative(rows)["assembly_accession"] == "GCA_903813345.2"


# ---------------------------------------------------------------- mitogenome provenance

def mrec(acc, bs="", length=15000, title="mitochondrion, complete genome", topo="circular", date="2020/01/01"):
    return {"accession": acc, "biosample": bs, "length": length, "title": title, "topology": topo,
            "completeness": "", "create_date": date}


def test_mito_in_assembly_wins():
    out = core.choose_mito("SAMA", set(), [("OX439053.1", 15641)], [mrec("OX439053.1", "SAMA"), mrec("MK1.1")])
    assert (out["mito_accession"], out["mito_source_derived"]) == ("OX439053.1", "same sample: in assembly")
    assert out["mito_substitute_candidates"] == "MK1.1"


def test_mito_same_biosample_then_linked_run_biosample():
    recs = [mrec("A1.1", "SAMX"), mrec("B1.1", "SAMA", length=9000), mrec("C1.1", "SAMH")]
    out = core.choose_mito("SAMA", {"SAMH"}, [], recs)
    assert out["mito_accession"] == "B1.1" and out["mito_source_derived"] == "same sample: same BioSample"
    assert "partial" in out["mito_note_derived"]           # same sample beats length; partial is noted
    out = core.choose_mito("SAMZ", {"SAMH"}, [], recs)
    assert out["mito_accession"] == "C1.1" and out["mito_source_derived"].startswith("same project sample")


def test_mito_substitute_flagged_and_prefers_insdc_complete():
    recs = [mrec("NC_081955.1"), mrec("MK858187.1"), mrec("XX1.1", "SAMQ", title="cox1 partial", topo="linear")]
    out = core.choose_mito("SAMA", set(), [], recs)
    assert out["mito_accession"] == "MK858187.1"           # INSDC original, not the RefSeq copy
    assert out["mito_source_derived"].startswith("SUBSTITUTE")
    assert "identity not verified" in out["mito_note_derived"]


def test_mito_short_other_sample_never_substitutes():
    out = core.choose_mito("SAMA", set(), [], [mrec("XX1.1", "SAMQ", length=600)])
    assert out["mito_accession"] == "" and out["mito_source_derived"] == ""


def test_mito_detected_by_name_when_ncbi_labels_it_chromosome():
    # live GCA_965119645.1 (Owenia): 18.8 kb mitogenome labelled 'Chromosome', chr_name 'mitochondrial_genome'
    reps = [{"genbank_accession": "OZ221967.1", "role": "assembled-molecule", "length": 18801,
             "assigned_molecule_location_type": "Chromosome", "chr_name": "mitochondrial_genome"}]
    assert core.mito_from_seqreport(reps) == [("OZ221967.1", 18801)]
    unnamed = [dict(reps[0], chr_name="1")]
    assert core.mito_from_seqreport(unnamed) == []
    assert core.mito_from_seqreport(unnamed, mito_only=True) == [("OZ221967.1", 18801)]


def test_mito_insdc_original_preferred_over_refseq_copy_even_if_copy_flagged_complete():
    copy = dict(mrec("NC_000931.1"), completeness="complete")
    orig = mrec("AF178678.1", title="Platynereis dumerilii mitochondrion, complete genome")
    assert core.choose_mito("SAMA", set(), [], [copy, orig])["mito_accession"] == "AF178678.1"


def test_unplaced_scaffold_named_mitochondrion_and_long_flag():
    # live GCA_965276055.1 (Genoscope): 'mitochondrion_1', role unplaced-scaffold, location 'Chromosome', 77 kb
    reps = [{"genbank_accession": "CBDHZN010000007.1", "role": "unplaced-scaffold", "length": 77074,
             "assigned_molecule_location_type": "Chromosome", "chr_name": "Un", "sequence_name": "mitochondrion_1"}]
    mt = core.mito_from_seqreport(reps)
    assert mt == [("CBDHZN010000007.1", 77074)]
    out = core.choose_mito("SAMA", set(), mt, [])
    assert out["mito_source_derived"] == "same sample: in assembly" and "verify" in out["mito_note_derived"]


def test_crosscheck_reports_differences():
    asm = [{"assembly_accession": "GCA_1.1", "paired_accession": "", "organism_name": "X", "assembly_level": "Scaffold",
            "biosample": "SAM1", "tax_id": "1", "found_in_derived": "Datasets+Entrez"},
           {"assembly_accession": "GCA_2.1", "paired_accession": "", "organism_name": "Y", "assembly_level": "Contig",
            "biosample": "SAM2", "tax_id": "2", "found_in_derived": "Datasets only"}]
    ent = [{"accession": "GCA_1.1", "level": "Chromosome", "biosample": "SAM1", "tax_id": "1"}]
    msgs = [q[3] for q in core.crosscheck(asm, ent)]
    assert any("assembly level: Datasets 'Scaffold' vs Entrez 'Chromosome'" in m for m in msgs)
    assert any("not in Entrez" in m for m in msgs)


def test_ftp_assembly_report_parsing_and_url():
    txt = ("# Assembly name:  ASM\n# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\t"
           "GenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
           "mitochondrion_1\tunplaced-scaffold\tna\tna\tCBDHZN010000007.1\t<>\tna\tPrimary Assembly\t77074\tna\n"
           "SUPER_1\tassembled-molecule\t1\tChromosome\tOZ000001.1\t<>\tna\tPrimary Assembly\t90000000\tna\n")
    seqs = core.parse_assembly_report(txt)
    assert core.mito_from_seqreport(seqs) == [("CBDHZN010000007.1", 77074)]
    assert core.ftp_report_url("GCA_949317995.1", "wpAchSqua1.1") == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/949/317/995/GCA_949317995.1_wpAchSqua1.1/"
        "GCA_949317995.1_wpAchSqua1.1_assembly_report.txt")
    assert core.ftp_report_url("x", "", "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/1/A_B").endswith("/A_B/A_B_assembly_report.txt")


def test_entrez_docsum_uses_genbank_accession_for_refseq_paired_record():
    j = {"result": {"uids": ["1"], "1": {"assemblyaccession": "GCF_000326865.1",
                                         "synonym": {"genbank": "GCA_000326865.1", "refseq": "GCF_000326865.1"}}}}
    assert core.entrez_assembly_docsums(j)[0]["accession"] == "GCA_000326865.1"


def test_incremental_cache_fetches_only_new_items(tmp_path):
    c = gi.Cache(tmp_path)
    calls = []

    def fetch(xs):
        calls.append(list(xs))
        return {x: x.upper() for x in xs}

    assert c.get_incremental("t", ["a", "b"], fetch, lambda a, b: {**a, **b}) == {"a": "A", "b": "B"}
    assert c.get_incremental("t", ["a", "b", "c"], fetch, lambda a, b: {**a, **b}) == {"a": "A", "b": "B", "c": "C"}
    assert c.get_incremental("t", ["a", "c"], fetch, lambda a, b: {**a, **b})["c"] == "C"
    assert calls == [["a", "b"], ["c"]]


# ---------------------------------------------------------------- same individual / confidence / GoaT / changes

def test_individual_evidence_acholoe_live_relations():
    # live ENA 2026-10-28: HiFi sample SAMEA13854475 and Hi-C sample SAMEA13854477 derived_from SAMEA13854387
    tsv = ("sample_accession\trelated_sample_accession\tscientific_name\n"
           "SAMEA13854475\tSAMEA13854387:derived_from\tAcholoe squamosa\n"
           "SAMEA13854387\tSAMEA13854475:has_part\tAcholoe squamosa\n"
           "SAMEA13854477\tSAMEA13854387:derived_from\tAcholoe squamosa\n")
    g = core.sample_graph(core.parse_ena_sample_relations(tsv))
    assert core.individual_evidence("SAMEA13854477", "SAMEA13854387", g, {}) == "same individual: BioSample relation"
    attrs = {"S1": {"tolid": "wpXxx1"}, "S2": {"tolid": "wpXxx2"}, "S3": {"tolid": "wpXxx1"}}
    assert core.individual_evidence("S1", "S2", {}, attrs).startswith("different individual")
    assert core.individual_evidence("S1", "S3", {}, attrs) == "same individual: same ToLID"
    assert core.individual_evidence("S9", "S2", {}, attrs) == "unknown"
    # two samples derived from a common specimen sample = same individual
    g2 = core.sample_graph({"A": [("P", "derived_from")], "B": [("P", "derived_from")]})
    assert core.individual_evidence("A", "B", g2, {}) == "same individual: BioSample relation"


def test_ena_sample_xml_tolid():
    xml = ("<SAMPLE_SET><SAMPLE accession='ERS1'><IDENTIFIERS><EXTERNAL_ID namespace='BioSample'>SAMEA13854387"
           "</EXTERNAL_ID></IDENTIFIERS><SAMPLE_ATTRIBUTES><SAMPLE_ATTRIBUTE><TAG>tolid</TAG><VALUE>wpAchSqua1</VALUE>"
           "</SAMPLE_ATTRIBUTE></SAMPLE_ATTRIBUTES></SAMPLE></SAMPLE_SET>")
    assert core.parse_ena_sample_xml(xml)["SAMEA13854387"]["tolid"] == "wpAchSqua1"


def test_category_confidence():
    pb = dict(Platform="PACBIO_SMRT")
    conf = lambda r, x=None: core.category_confidence(core.classify_run(r, x)[1])  # noqa: E731
    assert conf(run(Model="Sequel II", **pb), {"library_protocol": "PacBio - HiFi"}) == "high"
    assert conf(run(Model="Revio", **pb)) == "high"
    assert conf(run(Model="Sequel", **pb)) == "medium"
    assert conf(run(Model="Sequel IIe", **pb)) == "low"
    assert conf(run(Model="Sequel II", **pb)) == "low"


def test_change_report():
    old = [{"Assembly Accession": "GCA_1.1", "Organism Name": "X", "Assembly Level": "Scaffold"},
           {"Assembly Accession": "GCA_2.1", "Organism Name": "Y"}]
    new = [{"Assembly Accession": "GCA_1.1", "Organism Name": "X", "Assembly Level": "Chromosome"},
           {"Assembly Accession": "GCA_3.1", "Organism Name": "Z"}]
    ch = {(c["change"], c["assembly"], c.get("field", "")) for c in core.change_report(old, new)}
    assert ch == {("new assembly", "GCA_3.1", ""), ("no longer listed (suppressed / replaced?)", "GCA_2.1", ""),
                  ("field changed", "GCA_1.1", "Assembly Level")}


def test_goat_project_field_string_or_list():
    j = {"results": [{"result": {"taxon_id": 1, "scientific_name": "A", "lineage": [],
                                 "fields": {"long_list": {"value": "DTOL"}, "sequencing_status": {"value": "in_assembly"}}}},
                     {"result": {"taxon_id": 2, "scientific_name": "B", "lineage": [],
                                 "fields": {"long_list": {"value": ["ATLASEA", "DTOL"]}}}}]}
    assert [g["goat_projects"] for g in core.goat_species(j)] == ["DTOL", "ATLASEA;DTOL"]
