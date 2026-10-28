"""Pure functions (no network): parse NCBI responses, classify raw data, link runs to assemblies, build tables.

Derived columns end in `_derived`; the rule that produced each value is written next to it, so nothing is a
silent guess. Original NCBI values are never changed.
"""
import csv
import io
import json
import math
import re

from ncbi import xml_root

RANKS = ("phylum", "class", "subclass", "order", "suborder", "family", "genus", "species")
LEVEL_ORDER = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}

# Raw-data categories, in display order. Assembly-linked counts are reported per category.
CATEGORIES = ["HiFi", "CLR", "PacBio_unknown", "ONT", "Hi-C", "short_read_WGS", "10x_linked",
              "RNA-seq", "Iso-Seq", "ONT_RNA", "Other"]

HIC_RE = re.compile(r"hi-?c\b|omni-?c|arima|dovetail|phase genomics|proximo|proximity ligation|"
                    r"chromatin conformation|chromosome conformation", re.I)
HIFI_RE = re.compile(r"hifi|\bccs\b|circular consensus", re.I)
CLR_RE = re.compile(r"\bclr\b|continuous long read|subread", re.I)
TENX_RE = re.compile(r"10x|chromium|linked[- ]read|stlfr|tell-?seq", re.I)


# ---------------------------------------------------------------- parsing

def num(x):
    try:
        v = float(str(x).replace(",", ""))
        return 0.0 if math.isnan(v) else v
    except (TypeError, ValueError):
        return 0.0


def parse_runinfo(text):
    """SraRunInfo CSV -> [dict]; drops BOM, blank lines and repeated header rows (page concatenation)."""
    rows = []
    for r in csv.DictReader(io.StringIO(text.lstrip("﻿"))):
        run = (r.get("Run") or "").strip()
        if run and run != "Run":
            rows.append({k.lstrip("﻿"): (v or "") for k, v in r.items() if k is not None})
    return rows


def _t(el, path):
    x = el.find(path)
    return (x.text or "").strip() if x is not None and x.text else ""


def parse_sra_xml(xml):
    """efetch db=sra XML -> {run: {experiment_title, design_description, library_protocol, study_title,
    sample attributes of interest}}. Only text fields used for classification / provenance."""
    out = {}
    for pkg in xml_root(xml).iter("EXPERIMENT_PACKAGE"):
        d = {
            "experiment_title": _t(pkg, "EXPERIMENT/TITLE"),
            "design_description": _t(pkg, "EXPERIMENT/DESIGN/DESIGN_DESCRIPTION"),
            "library_protocol": _t(pkg, "EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_CONSTRUCTION_PROTOCOL"),
            "study_title": _t(pkg, "STUDY/DESCRIPTOR/STUDY_TITLE"),
        }
        for a in pkg.findall("SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
            tag = _t(a, "TAG").lower()
            if tag in ("tissue", "organism part", "specimen_voucher", "isolate", "sex", "geo_loc_name",
                       "geographic location (country and/or sea)", "collection_date", "collection date"):
                key = {"organism part": "tissue", "geographic location (country and/or sea)": "geo_loc_name",
                       "collection date": "collection_date"}.get(tag, tag)
                d.setdefault(key, _t(a, "VALUE"))
        for run in pkg.findall("RUN_SET/RUN"):
            out[run.get("accession", "")] = d
    return out


def parse_bioproject_parents(xml):
    """efetch db=bioproject XML -> {accession: [parent accessions]} (Hierarchical links, TopAdmin/TopSingle)."""
    out = {}
    for ds in xml_root(xml).iter("DocumentSummary"):
        acc_el = ds.find("Project/ProjectID/ArchiveID")
        acc = acc_el.get("accession", "") if acc_el is not None else ""
        parents = []
        for link in ds.findall("ProjectLinks/Link"):
            ref = link.find("ProjectIDRef")
            hier = link.find("Hierarchical")
            if ref is None or hier is None:
                continue
            if ref.get("accession") == acc:          # this project is the member -> MemberID is its parent
                for m in hier.findall("MemberID"):
                    if m.get("accession") and m.get("accession") != acc:
                        parents.append(m.get("accession"))
        if acc:
            out[acc] = sorted(set(parents))
    return out


def parse_taxonomy_xml(xml):
    """efetch db=taxonomy XML -> {taxid: {rank_derived, lineage, phylum..species}}; merged old taxids included."""
    out = {}
    for t in xml_root(xml).findall("Taxon"):
        r = {"rank_derived": t.findtext("Rank", ""), "lineage": t.findtext("Lineage", "")}
        r.update({k: "" for k in RANKS})
        ancestors = [x.findtext("TaxId", "") for x in t.findall("LineageEx/Taxon")]
        for x in t.findall("LineageEx/Taxon") + [t]:
            rank = x.findtext("Rank", "")
            if rank in RANKS:
                r[rank] = x.findtext("ScientificName", "")
        r["ancestors"] = ancestors
        tid = t.findtext("TaxId", "")
        out[tid] = r
        for aka in t.findall("AkaTaxIds/TaxId"):
            out[aka.text] = r
    return out


# Same column names and order as the NCBI Datasets web table / `dataformat tsv genome` export the user
# used before (Annelida_20260208.tsv), so old sheets and scripts keep working.
NCBI_COLS = [
    "Assembly Name", "Assembly Accession", "Assembly Paired Assembly Accession", "Organism Name",
    "Organism Taxonomic ID", "ANI Check status", "Organism Infraspecific Names Breed",
    "Organism Infraspecific Names Strain", "Organism Infraspecific Names Cultivar",
    "Organism Infraspecific Names Ecotype", "Organism Infraspecific Names Isolate",
    "Organism Infraspecific Names Sex", "Annotation Name", "Annotation Release Date",
    "Assembly Stats Total Sequence Length", "Assembly Stats Total Number of Chromosomes",
    "Assembly Stats Number of Contigs", "Assembly Level", "Assembly Release Date", "WGS project accession",
    "Assembly Stats Contig N50", "Assembly Stats Scaffold N50", "Assembly Stats Number of Scaffolds",
    "Assembly Stats GC Percent", "Annotation BUSCO Complete", "Annotation BUSCO Single Copy",
    "Annotation BUSCO Duplicated", "Annotation BUSCO Fragmented", "Annotation BUSCO Missing",
    "Annotation BUSCO Lineage", "Assembly Sequencing Tech", "Assembly Submitter",
    "Assembly BioProject Accession", "Assembly BioSample Accession", "Annotation Count Gene Total",
    "Annotation Count Gene Protein-coding", "Annotation Count Gene Pseudogene", "Type Material Display Text",
    "CheckM marker set", "CheckM completeness", "CheckM contamination"]
ANNOTATION_COLS = [c for c in NCBI_COLS if c.startswith("Annotation ")]


def _int(x):
    """'1481470071' -> 1481470071 so Excel sorts numbers as numbers; other values unchanged."""
    if isinstance(x, str) and x.isdigit():
        return int(x)
    return "" if x is None else x


def ncbi_columns(rep):
    """Datasets v2 genome report -> the NCBI web-table columns (NCBI_COLS), original values."""
    ai, st, org = rep.get("assembly_info", {}), rep.get("assembly_stats", {}), rep.get("organism", {})
    ann = rep.get("annotation_info") or {}
    busco, genes = ann.get("busco") or {}, (ann.get("stats") or {}).get("gene_counts") or {}
    inf, ck = org.get("infraspecific_names") or {}, rep.get("checkm_info") or {}
    v = [ai.get("assembly_name"), rep.get("accession"), rep.get("paired_accession"), org.get("organism_name"),
         org.get("tax_id"), (rep.get("average_nucleotide_identity") or {}).get("taxonomy_check_status"),
         inf.get("breed"), inf.get("strain"), inf.get("cultivar"), inf.get("ecotype"), inf.get("isolate"),
         inf.get("sex"), ann.get("name"), ann.get("release_date"),
         st.get("total_sequence_length"), st.get("total_number_of_chromosomes"), st.get("number_of_contigs"),
         ai.get("assembly_level"), ai.get("release_date"), (rep.get("wgs_info") or {}).get("wgs_project_accession"),
         st.get("contig_n50"), st.get("scaffold_n50"), st.get("number_of_scaffolds"), st.get("gc_percent"),
         busco.get("complete"), busco.get("single_copy"), busco.get("duplicated"), busco.get("fragmented"),
         busco.get("missing"), busco.get("busco_lineage"), ai.get("sequencing_tech"), ai.get("submitter"),
         ai.get("bioproject_accession"), (ai.get("biosample") or {}).get("accession"),
         genes.get("total"), genes.get("protein_coding"), genes.get("pseudogene"),
         (rep.get("type_material") or {}).get("type_display_text"),
         ck.get("checkm_marker_set"), ck.get("completeness"), ck.get("contamination")]
    return {c: _int(x) for c, x in zip(NCBI_COLS, v)}


def assembly_row(rep):
    """One Datasets v2 genome report -> flat assembly record: NCBI web-table columns + derived helpers."""
    ai, st = rep.get("assembly_info", {}), rep.get("assembly_stats", {})
    org, bs = rep.get("organism", {}), ai.get("biosample", {}) or {}
    attrs = {a.get("name"): a.get("value") for a in bs.get("attributes", []) if a.get("name")}
    ann = rep.get("annotation_info") or {}
    busco = ann.get("busco") or {}
    lineage = ai.get("bioproject_lineage") or []
    projects = {}
    for lin in lineage:
        for p in lin.get("bioprojects", []):
            projects[p["accession"]] = p.get("parent_accessions", [])
    bp = ai.get("bioproject_accession", "")
    infra = org.get("infraspecific_names") or {}
    return {
        **ncbi_columns(rep),
        # some "Complete Genome" assemblies are a single ~16 kb mitochondrion (e.g. GCA_965119645.1 Owenia)
        "assembly_content_derived": ("mitogenome-only (< 1 Mb)" if 0 < num(st.get("total_sequence_length")) < 1e6
                                     else "nuclear genome"),
        "organism_name": org.get("organism_name", ""),
        "tax_id": str(org.get("tax_id", "")),
        "assembly_accession": rep.get("accession", ""),
        "paired_accession": rep.get("paired_accession", ""),
        "assembly_name": ai.get("assembly_name", ""),
        "assembly_level": ai.get("assembly_level", ""),
        "assembly_type": ai.get("assembly_type", ""),
        "assembly_status": ai.get("assembly_status", ""),
        "refseq_category": ai.get("refseq_category", ""),
        "release_date": ai.get("release_date", ""),
        "submitter": ai.get("submitter", ""),
        "bioproject": bp,
        "bioproject_parents": ";".join(projects.get(bp, [])),
        "biosample": bs.get("accession", ""),
        "isolate": infra.get("isolate", "") or attrs.get("isolate", ""),
        "sex": infra.get("sex", "") or attrs.get("sex", ""),
        "tissue": attrs.get("tissue", "") or attrs.get("organism part", ""),
        "geo_loc_name": attrs.get("geo_loc_name", "") or attrs.get("geographic location (country and/or sea)", ""),
        "sequencing_tech": ai.get("sequencing_tech", ""),
        "assembly_method": ai.get("assembly_method", ""),
        "genome_coverage": st.get("genome_coverage", ""),
        "total_length_Mb": round(num(st.get("total_sequence_length")) / 1e6, 2) or "",
        "contig_N50_kb": round(num(st.get("contig_n50")) / 1e3, 1) or "",
        "scaffold_N50_kb": round(num(st.get("scaffold_n50")) / 1e3, 1) or "",
        "n_chromosomes": st.get("total_number_of_chromosomes", ""),
        "n_scaffolds": st.get("number_of_scaffolds", ""),
        "n_contigs": st.get("number_of_contigs", ""),
        "gc_percent": st.get("gc_percent", ""),
        "wgs_project": (rep.get("wgs_info") or {}).get("wgs_project_accession", ""),
        "annotation_provider": ann.get("provider", ""),
        "annotation_genes_total": (ann.get("stats") or {}).get("gene_counts", {}).get("total", ""),
        "annotation_BUSCO_complete": busco.get("complete", ""),
        "annotation_BUSCO_lineage": busco.get("busco_lineage", ""),
        "organelles_declared": ";".join(f"{o.get('description', '')}:{o.get('total_seq_length', '')}"
                                        for o in rep.get("organelle_info") or []),
        "linked_assemblies": ";".join(f"{x.get('linked_assembly', '')}({x.get('assembly_type', '')})"
                                      for x in ai.get("linked_assemblies") or []),
    }


def merge_gca_gcf(reports):
    """Keep one row per assembly: GCA rows, with the paired GCF in `paired_accession`.
    A GCF whose GCA is not in the result is kept as its own row (rare)."""
    rows = [assembly_row(r) for r in reports]
    gca = {r["assembly_accession"]: r for r in rows if r["assembly_accession"].startswith("GCA_")}
    out = list(gca.values())
    for r in rows:
        acc = r["assembly_accession"]
        if acc.startswith("GCF_"):
            pair = r["paired_accession"]
            if pair in gca:
                gca[pair]["paired_accession"] = acc
                # NCBI eukaryotic annotation is attached to the GCF; copy it onto the GCA row
                for k in ANNOTATION_COLS + ["annotation_provider", "annotation_genes_total",
                                            "annotation_BUSCO_complete", "annotation_BUSCO_lineage"]:
                    gca[pair][k] = gca[pair][k] if gca[pair][k] not in ("", None) else r[k]
                gca[pair]["Assembly Paired Assembly Accession"] = acc
            else:
                out.append(r)
    return out


def mito_from_seqreport(reports, mito_only=False):
    """Datasets sequence_reports -> [(genbank_accession, length)] of mitochondrial molecules in the assembly.
    NCBI labels are unreliable: GCA_965119645.1 has its 18.8 kb mitogenome as location type 'Chromosome' with
    chr_name 'mitochondrial_genome', so sequence names are checked too; for a mitogenome-only assembly
    (mito_only=True) every sequence is taken."""
    def is_mt(r):
        if (r.get("assigned_molecule_location_type") or "").lower() == "mitochondrion":
            return True
        names = f"{r.get('chr_name', '')} {r.get('sequence_name', '')}".lower()
        return bool(re.search(r"mito|\bmt\b|\bchrm\b", names))
    hits = [r for r in reports if is_mt(r)] or (list(reports) if mito_only else [])
    return [(r.get("genbank_accession", "") or r.get("refseq_accession", ""), r.get("length", "")) for r in hits]


def parse_assembly_report(text):
    """NCBI FTP *_assembly_report.txt -> sequence dicts with the Datasets sequence_reports keys
    (chr_name, sequence_name, role, assigned_molecule_location_type, genbank_accession, length)."""
    head, out = None, []
    for line in text.splitlines():
        if line.startswith("# Sequence-Name"):
            head = line[2:].split("\t")
            continue
        if line.startswith("#") or not line.strip() or head is None:
            continue
        d = dict(zip(head, line.split("\t")))
        na = lambda v: "" if v in ("na", None) else v  # noqa: E731
        out.append({"sequence_name": na(d.get("Sequence-Name")), "role": na(d.get("Sequence-Role")),
                    "chr_name": na(d.get("Assigned-Molecule")),
                    "assigned_molecule_location_type": na(d.get("Assigned-Molecule-Location/Type")),
                    "genbank_accession": na(d.get("GenBank-Accn")), "length": _int(na(d.get("Sequence-Length")))})
    return out


def ftp_report_url(accession, assembly_name, ftp_path=""):
    """URL of the FTP assembly report; Entrez's ftppath_genbank when known, else built from the accession."""
    if ftp_path:
        base = ftp_path.replace("ftp://", "https://")
    else:
        n = accession[4:13]
        safe = re.sub(r"[^A-Za-z0-9._-]", "_", assembly_name.strip())
        base = (f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{accession[:3]}/{n[:3]}/{n[3:6]}/{n[6:9]}/"
                f"{accession}_{safe}")
    return f"{base}/{base.rstrip('/').rsplit('/', 1)[-1]}_assembly_report.txt"


def mito_docsums(esummary_json):
    """esummary db=nuccore JSON -> [record dict]."""
    res = esummary_json.get("result", {})
    out = []
    for uid in res.get("uids", []):
        x = res[uid]
        out.append({"accession": x.get("accessionversion") or x.get("caption", ""),
                    "title": x.get("title", ""), "length": x.get("slen", ""), "tax_id": str(x.get("taxid", "")),
                    "organism": x.get("organism", ""), "biosample": x.get("biosample", "") or "",
                    "topology": x.get("topology", ""), "completeness": x.get("completeness", "") or "",
                    "create_date": x.get("createdate", ""), "source_db": x.get("sourcedb", "")})
    return out


# ---------------------------------------------------------------- classification

def classify_run(r, x=None):
    """Raw-data category of one run. r = RunInfo row, x = SRA XML text fields (may be empty).
    Returns (category, rule). Keyword evidence is searched in experiment title, design description,
    library name and library construction protocol (never the study abstract)."""
    x = x or {}
    strat = r.get("LibraryStrategy", "").upper()
    src = r.get("LibrarySource", "").upper()
    plat = r.get("Platform", "").upper()
    model = r.get("Model", "")
    sel = r.get("LibrarySelection", "")
    text = " ".join([x.get("experiment_title", ""), x.get("design_description", ""),
                     x.get("library_protocol", ""), r.get("LibraryName", "")])
    rna = src in ("TRANSCRIPTOMIC", "TRANSCRIPTOMIC SINGLE CELL") or strat in ("RNA-SEQ", "MIRNA-SEQ", "SSRNA-SEQ",
                                                                             "NCRNA-SEQ", "FL-CDNA", "EST")
    if src in ("METAGENOMIC", "METATRANSCRIPTOMIC"):
        return "Other", f"LibrarySource={src}"
    if strat == "HI-C":
        return "Hi-C", "LibraryStrategy=Hi-C"
    if not rna and (HIC_RE.search(sel) or HIC_RE.search(text)):
        m = HIC_RE.search(sel) or HIC_RE.search(text)
        return "Hi-C", f"keyword '{m.group(0)}' ({strat or 'no strategy'})"
    if plat == "PACBIO_SMRT":
        if rna:
            return "Iso-Seq", f"PacBio & {src or strat}"
        if HIFI_RE.search(text):
            return "HiFi", f"PacBio & keyword '{HIFI_RE.search(text).group(0)}'"
        if CLR_RE.search(text):
            return "CLR", f"PacBio & keyword '{CLR_RE.search(text).group(0)}'"
        if re.search(r"revio|vega", model, re.I):
            return "HiFi", f"Model={model} (HiFi-only instrument)"
        if re.search(r"sequel iie", model, re.I):
            return "HiFi", f"Model={model}, no keyword (on-board CCS instrument; inferred)"
        if re.fullmatch(r"(pacbio )?rs( ii)?|sequel", model.strip(), re.I):
            return "CLR", f"Model={model} (pre-HiFi instrument)"
        return "PacBio_unknown", f"Model={model}, no HiFi/CLR keyword"
    if plat == "OXFORD_NANOPORE":
        return ("ONT_RNA", f"ONT & {src or strat}") if rna else ("ONT", "Platform=OXFORD_NANOPORE")
    if rna:
        return "RNA-seq", f"short read & {src or strat}"
    if strat in ("WGS", "WGA", "OTHER", "") and src in ("GENOMIC", "GENOMIC SINGLE CELL", ""):
        m = TENX_RE.search(text)
        if m:
            return "10x_linked", f"keyword '{m.group(0)}'"
        if strat == "WGS" and plat != "CAPILLARY":
            return "short_read_WGS", f"WGS & {plat}"
    return "Other", f"LibraryStrategy={r.get('LibraryStrategy', '')}"


# ---------------------------------------------------------------- linking

def link_evidence(run, asm, bp_parents):
    """How a run (same species) is connected to an assembly.
    biosample > bioproject (same project or parent/child) > sibling_project (share a direct parent) > ''."""
    if run.get("BioSample") and run["BioSample"] == asm["biosample"]:
        return "biosample"
    rbp, abp = run.get("BioProject", ""), asm["bioproject"]
    if not rbp or not abp:
        return ""
    a_par = set(filter(None, asm["bioproject_parents"].split(";"))) | set(bp_parents.get(abp, []))
    r_par = set(bp_parents.get(rbp, []))
    if rbp == abp or rbp in a_par or abp in r_par:
        return "bioproject"
    if a_par & r_par:
        return "sibling_project"
    return ""


def representative(asms):
    """Best assembly of one species: nuclear genome > primary (not alternate) > RefSeq reference > assembly level
    > contig N50 > newest."""
    def key(a):
        alt = "alternate" in a["assembly_type"].lower()
        return (a["assembly_content_derived"] == "nuclear genome", not alt, a["refseq_category"] == "reference genome", LEVEL_ORDER.get(a["assembly_level"], 0),
                num(a["contig_N50_kb"]), a["release_date"])
    return max(asms, key=key)


def declared_flags(tech):
    """Technologies the submitter declared in `sequencing_tech` -> set of categories (for QC only)."""
    t = tech or ""
    f = set()
    if HIFI_RE.search(t) or re.search(r"pacbio|sequel|revio", t, re.I):
        f.add("PacBio")
    if re.search(r"nanopore|\bont\b|promethion|minion|gridion", t, re.I):
        f.add("ONT")
    if HIC_RE.search(t):
        f.add("Hi-C")
    return f


def join(xs, limit=None):
    xs = [x for x in dict.fromkeys(xs) if x]
    if limit and len(xs) > limit:
        return ";".join(xs[:limit]) + f";...(+{len(xs) - limit})"
    return ";".join(xs)


def to_json(x):
    return json.dumps(x, ensure_ascii=False)


def norm_taxid(s):
    """'46593.0' / ' 46593 ' -> '46593'; empty -> ''."""
    s = str(s or "").strip()
    if s.lower() in ("", "nan", "none"):
        return ""
    return s[:-2] if s.endswith(".0") else s


# ---------------------------------------------------------------- mitogenome provenance

MITO_FULL = 10000      # shorter records are reported as partial
MITO_LONG = 30000      # annelid mitogenomes are ~15-23 kb; longer "mitochondrion" sequences are flagged
MITO_SOURCES = {
    "in_assembly": "same sample: in assembly",
    "same_biosample": "same sample: same BioSample",
    "linked_run_biosample": "same project sample: BioSample of linked raw data",
    "substitute": "SUBSTITUTE: same species name, different sample - verify identity",
}


def _base(acc):
    return str(acc).split(".")[0]


def _mito_key(m):
    """Best record within one provenance level: INSDC original (not a RefSeq NC_ copy) > complete > circular
    > length >= 10 kb > newest."""
    complete = m.get("completeness") == "complete" or "complete genome" in m.get("title", "").lower()
    return (not str(m.get("accession", "")).startswith("NC_"), complete, m.get("topology") == "circular",
            num(m.get("length")) >= MITO_FULL, m.get("create_date", ""))


def choose_mito(asm_biosample, linked_biosamples, in_assembly, records):
    """One mitogenome per assembly, from the same sample whenever possible.

    in_assembly: [(accession, length)] mitochondrial molecules of the assembly itself.
    records: nuccore mitochondrion records of the species plus those found by BioSample (dicts).
    Order: in_assembly > same_biosample > linked_run_biosample > substitute (same species name, different
    sample: NCBI identifications are not checked, so these are flagged for verification).
    """
    by_base = {_base(m["accession"]): m for m in records}
    inside = {_base(a) for a, _ in in_assembly}
    levels = {"in_assembly": [dict(by_base.get(_base(a), {}), accession=a, length=n) for a, n in in_assembly]}
    rest = [m for m in records if _base(m["accession"]) not in inside]
    levels["same_biosample"] = [m for m in rest if asm_biosample and m.get("biosample") == asm_biosample]
    other_bs = set(linked_biosamples) - {asm_biosample, ""}
    levels["linked_run_biosample"] = [m for m in rest if m.get("biosample") in other_bs]
    same = {_base(m["accession"]) for k in ("same_biosample", "linked_run_biosample") for m in levels[k]}
    levels["substitute"] = [m for m in rest if _base(m["accession"]) not in same and num(m.get("length")) >= MITO_FULL]

    out = {"mito_accession": "", "mito_length_bp": "", "mito_source_derived": "", "mito_note_derived": "",
           "mito_same_sample_all": join(m["accession"] for k in ("in_assembly", "same_biosample",
                                                                   "linked_run_biosample") for m in levels[k]),
           "mito_substitute_n": len(levels["substitute"]) or "",
           "mito_substitute_candidates": join((m["accession"] for m in
                                               sorted(levels["substitute"], key=_mito_key, reverse=True)), 20),
           "_levels": {k: [m["accession"] for m in v] for k, v in levels.items()}}
    for k in MITO_SOURCES:
        if levels[k]:
            best = max(levels[k], key=_mito_key) if k != "in_assembly" else levels[k][0]
            notes = []
            if num(best.get("length")) < MITO_FULL:
                notes.append(f"partial (< {MITO_FULL} bp)")
            if num(best.get("length")) > MITO_LONG:
                notes.append(f"> {MITO_LONG} bp: possibly concatenated copies or mis-labelled, verify")
            if k == "substitute":
                bs = {m.get("biosample") or m["accession"] for m in levels[k]}
                notes.append(f"{len(levels[k])} same-name record(s) from {len(bs)} sample(s); "
                             "identity not verified")
            out.update(mito_accession=best["accession"], mito_length_bp=_int(best.get("length", "")),
                       mito_source_derived=MITO_SOURCES[k], mito_note_derived="; ".join(notes))
            break
    return out


# ---------------------------------------------------------------- second source: Entrez assembly db

def entrez_assembly_docsums(esummary_json):
    res = esummary_json.get("result", {})
    out = []
    for uid in res.get("uids", []):
        x = res[uid]
        syn = x.get("synonym") or {}
        # for a RefSeq-paired record assemblyaccession can be the GCF (live: Helobdella GCF_000326865.1);
        # the GenBank accession is in synonym.genbank
        out.append({"accession": syn.get("genbank") or x.get("assemblyaccession", ""), "refseq": syn.get("refseq", ""),
                    "name": x.get("assemblyname", ""), "species": x.get("speciesname", ""),
                    "organism": x.get("organism", ""), "tax_id": str(x.get("taxid", "")),
                    "level": x.get("assemblystatus", ""), "biosample": x.get("biosampleaccn", ""),
                    "bioproject": ";".join(p.get("bioprojectaccn", "") for p in x.get("gb_bioprojects") or []),
                    "submitter": x.get("submitterorganization", ""),
                    "release": (x.get("asmreleasedate_genbank") or "")[:10].replace("/", "-"),
                    "ftp": x.get("ftppath_genbank", "")})
    return out


def crosscheck(asm_rows, entrez):
    """QC rows comparing Datasets assemblies with Entrez assembly docsums."""
    qc = []
    ent = {e["accession"]: e for e in entrez}
    for a in asm_rows:
        acc = a["assembly_accession"]
        e = ent.get(acc) or ent.get(a["paired_accession"])
        if a["found_in_derived"].startswith("Entrez only"):
            qc.append(("crosscheck", acc, a["organism_name"], "missing from the Datasets taxon query, present in "
                       "Entrez assembly db; Datasets report fetched by accession"))
        if e is None:
            if acc.startswith("GCA_"):
                qc.append(("crosscheck", acc, a["organism_name"], "in Datasets but not in Entrez assembly "
                           "(latest_genbank) - may be suppressed/replaced or not yet indexed"))
            continue
        for field, dv, ev in (("assembly level", a["assembly_level"], e["level"]),
                              ("BioSample", a["biosample"], e["biosample"]),
                              ("taxid", a["tax_id"], e["tax_id"])):
            if dv and ev and str(dv) != str(ev):
                qc.append(("crosscheck", acc, a["organism_name"], f"{field}: Datasets '{dv}' vs Entrez '{ev}'"))
    return qc


# ---------------------------------------------------------------- same individual? (ENA sample relations / ToLID)

def parse_ena_sample_relations(tsv):
    """ENA portal result=sample (sample_accession, related_sample_accession) -> {acc: [(related, type)]}.
    related_sample_accession looks like 'SAMEA13854387:derived_from;SAMEA1:same_as'."""
    out = {}
    lines = [x for x in tsv.splitlines() if x.strip()]
    if not lines:
        return out
    head = lines[0].split("\t")
    for line in lines[1:]:
        d = dict(zip(head, line.split("\t")))
        rel = []
        for item in filter(None, (d.get("related_sample_accession") or "").split(";")):
            acc, _, typ = item.strip().partition(":")
            if acc:
                rel.append((acc, typ))
        out[d.get("sample_accession", "")] = rel
    return out


def parse_ena_sample_xml(xml):
    """ENA browser XML of samples -> {biosample: {"tolid": .., "specimen_id": ..}}."""
    out = {}
    for s in xml_root(xml).iter("SAMPLE"):
        bs = ""
        for e in s.findall("IDENTIFIERS/EXTERNAL_ID"):
            if e.get("namespace") == "BioSample":
                bs = (e.text or "").strip()
        bs = bs or s.get("accession", "")
        attrs = {(_t(a, "TAG")).lower(): _t(a, "VALUE") for a in s.findall("SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")}
        out[bs] = {"tolid": attrs.get("tolid", ""), "specimen_id": attrs.get("specimen_id", "")}
    return out


def sample_graph(relations):
    """{acc: [(related, type)]} -> undirected adjacency {acc: set(accs)} (relations are recorded on one side)."""
    g = {}
    for a, rel in relations.items():
        for b, _ in rel:
            g.setdefault(a, set()).add(b)
            g.setdefault(b, set()).add(a)
    return g


def individual_evidence(run_bs, asm_bs, graph, attrs):
    """Is the run's sample the assembly's individual? Returns 'same individual: <basis>',
    'different individual: <basis>' or 'unknown'. graph = sample_graph(...); attrs = {biosample: {tolid}}.
    A ToLID names one individual (wpAchSqua1 vs wpAchSqua2)."""
    if not run_bs or not asm_bs:
        return "unknown"
    if run_bs == asm_bs:
        return "same individual: same BioSample"
    n_run, n_asm = graph.get(run_bs, set()), graph.get(asm_bs, set())
    if asm_bs in n_run or (n_run & n_asm):
        return "same individual: BioSample relation"
    t1, t2 = (attrs.get(run_bs) or {}).get("tolid", ""), (attrs.get(asm_bs) or {}).get("tolid", "")
    if t1 and t2:
        return "same individual: same ToLID" if t1 == t2 else f"different individual: ToLID {t1} vs {t2}"
    return "unknown"


def category_confidence(rule):
    """Confidence of category_derived from its rule text: high = explicit strategy/keyword or HiFi-only
    instrument; medium = pre-HiFi instrument; low = inferred or unresolved."""
    r = rule.lower()
    if "inferred" in r or "no hifi/clr keyword" in r:
        return "low"
    if "pre-hifi instrument" in r:
        return "medium"
    return "high"


# ---------------------------------------------------------------- GoaT (third source) / change report

def goat_assemblies(j):
    return [{"accession": x["result"].get("assembly_id", ""), "scientific_name": x["result"].get("scientific_name", ""),
             "tax_id": str(x["result"].get("taxon_id", ""))} for x in j.get("results", [])]


def _as_list(v):
    """GoaT list fields are a list for several values but a plain string for one."""
    if not v:
        return ""
    return ";".join(v) if isinstance(v, list) else str(v)


def goat_species(j):
    out = []
    for x in j.get("results", []):
        r = x["result"]
        lin = {n.get("taxon_rank"): n.get("scientific_name") for n in r.get("lineage", [])}
        out.append({"tax_id": str(r.get("taxon_id", "")), "scientific_name": r.get("scientific_name", ""),
                    "class": lin.get("class", ""), "order": lin.get("order", ""), "family": lin.get("family", ""),
                    "goat_sequencing_status": ((r.get("fields") or {}).get("sequencing_status") or {}).get("value", ""),
                    # projects that list the species (DTOL, ERGA, ATLASEA, ...): GoaT field long_list
                    "goat_projects": _as_list(((r.get("fields") or {}).get("long_list") or {}).get("value"))})
    return out


CHANGE_FIELDS = ["Organism Name", "Assembly Level", "Assembly Stats Contig N50", "Assembly Stats Scaffold N50",
                 "Assembly BioSample Accession", "data_types_derived", "mito_accession", "mito_source_derived"]


def change_report(old_rows, new_rows):
    """Compare the previous run's Assemblies table with the new one."""
    old = {r["Assembly Accession"]: r for r in old_rows}
    new = {str(r["Assembly Accession"]): r for r in new_rows}
    out = []
    for acc in sorted(set(new) - set(old)):
        out.append({"change": "new assembly", "assembly": acc, "organism": new[acc]["Organism Name"]})
    for acc in sorted(set(old) - set(new)):
        out.append({"change": "no longer listed (suppressed / replaced?)", "assembly": acc,
                    "organism": old[acc].get("Organism Name", "")})
    for acc in sorted(set(old) & set(new)):
        for f in CHANGE_FIELDS:
            a, b = str(old[acc].get(f, "") or ""), str(new[acc].get(f, "") or "")
            if a != b:
                out.append({"change": "field changed", "assembly": acc, "organism": new[acc]["Organism Name"],
                            "field": f, "old": a, "new": b})
    return out


def tsv_safe(v):
    """Neutralise spreadsheet formula injection in TSV cells (submitter text such as titles goes into the
    files): a text cell starting with = + - @ TAB or CR gets a leading apostrophe. Numbers, including
    negative ones, are left as they are. The .xlsx output writes inline strings and needs no escaping."""
    if not isinstance(v, str) or not v or v[0] not in "=+-@\t\r":
        return v
    try:
        float(v)
        return v
    except ValueError:
        return "'" + v
