#!/usr/bin/env python3
"""One-command inventory of NCBI whole-genome assemblies for a taxon.

For every genome assembly under a taxon (e.g. Annelida) it collects:
  - assembly accession / level / N50 / submitter / BioProject / BioSample (NCBI Datasets v2)
  - the raw data behind it: HiFi, CLR, ONT, Hi-C, short-read WGS, 10x, RNA-seq, Iso-Seq ... (SRA RunInfo + XML),
    linked by BioSample, BioProject, or sibling BioProject under the same umbrella project
  - one mitogenome per assembly from the same sample (inside the assembly > same BioSample > BioSample of the
    linked raw data); only if none exists a same-name record is used and flagged SUBSTITUTE
  - taxonomy: NCBI lineage (phylum ... species) and, optionally, the WoRMS accepted name and classification

Usage:
  python3 genome_inventory.py Annelida -o Annelida_genomes
  python3 genome_inventory.py Polynoidae --no-worms
  python3 genome_inventory.py 6340 -o out            # taxid also accepted
Raw responses are cached in <out>/raw/; rerunning resumes from the cache, --refresh refetches everything.
Standard library only (Python >= 3.8). Optional: NCBI_EMAIL / NCBI_API_KEY environment variables.
"""
import argparse
import csv
import datetime
import json
import math
import sys
import urllib.parse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import core  # noqa: E402
from ncbi import Client, NcbiError, resolve_taxid  # noqa: E402
from xlsx import write_xlsx  # noqa: E402

TAXID_BATCH = 100      # OR clauses per Entrez query
MITO_MIN, MITO_MAX = 10000, 100000   # nuccore length window treated as (near-)complete mitogenome


def log(msg):
    print(f"[{datetime.datetime.now():%H:%M:%S}] {msg}", flush=True)


class Cache:
    def __init__(self, d, refresh=False):
        self.d = Path(d)
        self.d.mkdir(parents=True, exist_ok=True)
        self.refresh = refresh

    def get_incremental(self, name, items, fn, merge):
        """Cache keyed by the input items (taxids, BioSamples ...): when NCBI adds new assemblies/species,
        only the new items are fetched and merged, so an update never serves stale results."""
        f = self.d / f"{name}.json"
        items = sorted(set(items))
        old = json.loads(f.read_text(encoding="utf-8")) if f.exists() and not self.refresh else None
        if old is not None and not (isinstance(old, dict) and "items" in old and "data" in old):
            old = None                                    # cache written by an older version: refetch
        if old is None:
            data, done = fn(items), set()
        else:
            done, data = set(old["items"]), old["data"]
            new = [x for x in items if x not in done]
            if not new:
                return data
            log(f"  {name}: {len(new)} new item(s) since the cached run -> fetching them")
            data = merge(data, fn(new))
        tmp = f.with_suffix(".tmp")
        tmp.write_text(json.dumps({"items": sorted(done | set(items)), "data": data}, ensure_ascii=False),
                       encoding="utf-8")
        tmp.replace(f)
        return data

    def get(self, name, fn, refresh=False):
        f = self.d / f"{name}.json"
        if f.exists() and not (self.refresh or refresh):
            return json.loads(f.read_text(encoding="utf-8"))
        v = fn()
        tmp = f.with_suffix(".tmp")
        tmp.write_text(json.dumps(v, ensure_ascii=False), encoding="utf-8")
        tmp.replace(f)
        return v


def batches(xs, n):
    xs = list(xs)
    for i in range(0, len(xs), n):
        yield xs[i:i + n]


# ---------------------------------------------------------------- fetch steps

def fetch_assemblies(cl, taxid, all_versions):
    params = {} if all_versions else {"filters.assembly_version": "current"}
    reps = list(cl.datasets_pages(f"genome/taxon/{taxid}/dataset_report", params))
    log(f"Datasets: {len(reps)} assembly reports")
    return reps


def fetch_entrez_assemblies(cl, taxid, all_versions):
    """Second, independent source: Entrez `assembly` db (NCBI keeps E-utilities running after the Assembly
    website was retired in 2024). Used to cross-check the Datasets list."""
    term = f"txid{taxid}[Organism:exp]" + ("" if all_versions else " AND latest_genbank[filter]")
    h = cl.esearch("assembly", term)
    out = []
    for j in cl.esummary_history("assembly", h):
        out += core.entrez_assembly_docsums(j)
    log(f"Entrez assembly: {len(out)} records")
    return out


def fetch_assemblies_by_accession(cl, accs):
    reps = []
    for b in batches(sorted(accs), 50):
        reps += list(cl.datasets_pages(f"genome/accession/{','.join(b)}/dataset_report"))
    return reps


SEQREPORT_MAX_PAGES = 20    # Datasets fallback only
SCAN_WORKERS = 6            # parallel FTP downloads (FTP files are not under the E-utilities rate limit)


def scan_assembly_mito(cl, rep, ftp_path, mito_only):
    """Mitochondrial molecules of one assembly. Method 1: FTP assembly_report (one request, all sequences,
    names included - e.g. Genoscope GCA_965276055.1 'mitochondrion_1' as an unplaced scaffold).
    Method 2 (fallback): Datasets sequence_reports. Returns (molecules, method or error)."""
    acc = rep["accession"]
    errors = []
    try:
        txt = cl.get_text(core.ftp_report_url(acc, rep["assembly_info"].get("assembly_name", ""), ftp_path))
        seqs = core.parse_assembly_report(txt)
        if seqs:
            return core.mito_from_seqreport(seqs, mito_only), "ftp assembly_report"
        errors.append("ftp: empty report")
    except NcbiError as e:
        errors.append(f"ftp: {e}"[:150])
    try:
        seqs = list(cl.datasets_pages(f"genome/accession/{acc}/sequence_reports", max_pages=SEQREPORT_MAX_PAGES))
        return core.mito_from_seqreport(seqs, mito_only), "datasets sequence_reports"
    except NcbiError as e:
        errors.append(f"datasets: {e}"[:150])
    return [], "ERROR " + " | ".join(errors)


def fetch_mito_in_assembly(cl, reps, ftp, mito_only, cache_file):
    """Scan every assembly; progress saved every 20 assemblies so a rerun continues where it stopped."""
    cache_file = Path(cache_file)
    state = json.loads(cache_file.read_text(encoding="utf-8")) if cache_file.exists() else {}
    if state and not isinstance(next(iter(state.values())), dict):        # cache of an older version
        state = {}
    todo = [r for r in reps if r["accession"] not in state or state[r["accession"]]["method"].startswith("ERROR")]
    with ThreadPoolExecutor(max_workers=SCAN_WORKERS) as ex:
        futs = {ex.submit(scan_assembly_mito, cl, r, ftp.get(r["accession"], ""), r["accession"] in mito_only):
                r["accession"] for r in todo}
        for i, fut in enumerate(as_completed(futs), 1):
            mol, method = fut.result()
            state[futs[fut]] = {"mito": mol, "method": method}
            if i % 20 == 0 or i == len(todo):
                cache_file.write_text(json.dumps(state), encoding="utf-8")
                log(f"  assembly sequence lists {i}/{len(todo)}")
    return state


def fetch_taxonomy(cl, taxids):
    out = {}
    for b in batches(sorted(taxids), 200):
        out.update(core.parse_taxonomy_xml(cl.efetch_ids("taxonomy", b)))
    return out


def fetch_sra(cl, taxids):
    """RunInfo + XML text fields for all SRA runs of the given species (incl. lower ranks)."""
    rows, fields = [], {}
    for b in batches(sorted(taxids), TAXID_BATCH):
        term = " OR ".join(f"txid{t}[Organism:exp]" for t in b)
        h = cl.esearch("sra", term)
        log(f"  SRA: {h['count']} experiments for {len(b)} taxa")
        if not h["count"]:
            continue
        for page in cl.efetch_history("sra", h, rettype="runinfo", retmode="text"):
            rows += core.parse_runinfo(page)
        for page in cl.efetch_history("sra", h, retmode="xml", page=200):
            fields.update(core.parse_sra_xml(page))
    uniq = {r["Run"]: r for r in rows}
    return {"runinfo": list(uniq.values()), "xml": fields}


def fetch_bioprojects(cl, accs):
    out = {}
    for b in batches(sorted(accs), TAXID_BATCH):
        h = cl.esearch("bioproject", " OR ".join(f"{a}[Project Accession]" for a in b))
        for page in cl.efetch_history("bioproject", h, retmode="xml", page=200):
            out.update(core.parse_bioproject_parents(page))
    return out


def fetch_mitogenomes(cl, taxids):
    recs = []
    for b in batches(sorted(taxids), TAXID_BATCH):
        term = (f"({' OR '.join(f'txid{t}[Organism:exp]' for t in b)}) AND mitochondrion[filter] "
                f"AND {MITO_MIN}:{MITO_MAX}[SLEN]")
        h = cl.esearch("nuccore", term)
        for j in cl.esummary_history("nuccore", h):
            recs += core.mito_docsums(j)
    return recs


def fetch_mito_by_biosample(cl, biosamples):
    """Mitochondrial nuccore records of the given BioSamples, any length (partial same-sample records count)."""
    recs = []
    for b in batches(sorted(biosamples), TAXID_BATCH):
        h = cl.esearch("nuccore", f"({' OR '.join(f'{x}[BioSample]' for x in b)}) AND mitochondrion[filter]")
        for j in cl.esummary_history("nuccore", h):
            recs += core.mito_docsums(j)
    return recs


ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/search?"
ENA_XML = "https://www.ebi.ac.uk/ena/browser/api/xml/"
GOAT = "https://goat.genomehubs.org/api/v2/search?"


def fetch_sample_info(cl, biosamples):
    """ENA sample relations (derived_from / same_as / has_part) and ToLID of each BioSample. A failed batch is
    skipped (logged), never fatal: the evidence then stays 'unknown'."""
    rel, attrs = {}, {}
    for b in batches(sorted(biosamples), 40):
        q = urllib.parse.urlencode({"result": "sample", "format": "tsv",
                                    "fields": "sample_accession,related_sample_accession",
                                    "query": " OR ".join(f'sample_accession="{x}"' for x in b)})
        try:
            rel.update(core.parse_ena_sample_relations(cl.get_text(ENA_PORTAL + q)))
        except NcbiError as e:
            log(f"  ENA sample relations batch skipped: {e}"[:200])
    def xml_batch(b):
        # one accession unknown to ENA makes the whole request 404: split the batch and retry the halves
        try:
            attrs.update(core.parse_ena_sample_xml(cl.get_text(ENA_XML + ",".join(b))))
        except (NcbiError, ValueError) as e:
            if len(b) > 1:
                xml_batch(b[:len(b) // 2])
                xml_batch(b[len(b) // 2:])
            else:
                log(f"  ENA sample XML skipped for {b[0]}: {e}"[:160])

    # ToLIDs exist only on ENA/ToL samples (SAMEA); NCBI-origin SAMN samples often 404 at ENA
    for b in batches(sorted(x for x in biosamples if x.startswith("SAMEA")), 50):
        xml_batch(b)
    return {"relations": rel, "attrs": attrs}


def fetch_goat(cl, taxid):
    """GoaT (Genomes on a Tree): assemblies and per-species sequencing status (third, independent source)."""
    def q(**p):
        return json.loads(cl.get_text(GOAT + urllib.parse.urlencode(p).replace("+", "%20")))
    a = q(query=f"tax_tree({taxid})", result="assembly", size=10000, fields="assembly_level")
    t = q(query=f"tax_tree({taxid}) AND tax_rank(species)", result="taxon", size=10000, fields="sequencing_status,long_list")
    return {"assemblies": core.goat_assemblies(a), "species": core.goat_species(t)}


def fetch_worms(cl, names):
    out = {}
    for b in batches(sorted(names), 50):
        for name, cands in zip(b, cl.worms_match(b)):
            c = (cands or [None])[0]
            out[name] = c or {}
    return out


# ---------------------------------------------------------------- build tables

def species_of(run_taxid, asm_taxids, tax):
    """Assembly taxa a run belongs to: same taxid, or the assembly taxid is an ancestor (subspecies...)."""
    anc = set(tax.get(run_taxid, {}).get("ancestors", []))
    return [t for t in asm_taxids if t == run_taxid or t in anc]


def build(taxon, taxid, asm_reports, mito_asm, tax, sra, bp_parents, mitos, worms, args, entrez=(), added=()):
    asms = core.merge_gca_gcf(asm_reports)
    ent = {e["accession"]: e for e in entrez}
    for a in asms:
        acc = a["assembly_accession"]
        a["found_in_derived"] = ("Entrez only (Datasets taxon query missed it; report fetched by accession)"
                                 if acc in added else "Datasets+Entrez" if acc in ent or a["paired_accession"] in ent
                                 else "Datasets only")
    asm_taxids = sorted({a["tax_id"] for a in asms})
    by_tax = defaultdict(list)
    for a in asms:
        by_tax[a["tax_id"]].append(a)

    # --- runs: classify + link
    runs = []
    for r in sra["runinfo"]:
        x = sra["xml"].get(r["Run"], {})
        cat, rule = core.classify_run(r, x)
        for t in species_of(core.norm_taxid(r.get("TaxID")), asm_taxids, tax) or []:
            links = [(a["assembly_accession"], ev) for a in by_tax[t]
                     if (ev := core.link_evidence(r, a, bp_parents))]
            runs.append({
                "species_tax_id": t, "Run": r["Run"], "ScientificName": r.get("ScientificName", ""),
                "TaxID": r.get("TaxID", ""), "category_derived": cat, "category_rule_derived": rule,
                "linked_assemblies_derived": ";".join(f"{acc}:{ev}" for acc, ev in links),
                "bases_Gb": round(core.num(r.get("bases")) / 1e9, 3),
                "avgLength": r.get("avgLength", ""), "spots": r.get("spots", ""),
                **{k: r.get(k, "") for k in ("Experiment", "BioProject", "BioSample", "LibraryName",
                                             "LibraryStrategy", "LibrarySelection", "LibrarySource",
                                             "LibraryLayout", "Platform", "Model", "ReleaseDate", "CenterName")},
                **{k: x.get(k, "") for k in ("experiment_title", "library_protocol", "design_description",
                                             "study_title", "tissue", "specimen_voucher", "isolate", "sex",
                                             "geo_loc_name", "collection_date")},
            })

    # --- mitogenome records -> species (relation to assemblies is filled after choose_mito)
    mito_rows = []
    for m in mitos:
        mito_rows.append(dict(m, species_tax_id=";".join(species_of(m["tax_id"], asm_taxids, tax)),
                              relation_derived=""))

    # --- taxonomy columns
    def taxo(t):
        n = tax.get(t, {})
        d = {k: n.get(k, "") for k in core.RANKS}
        d["ncbi_rank"] = n.get("rank_derived", "")
        return d

    def worms_cols(name):
        w = worms.get(name) or {}
        return {"worms_valid_name": w.get("valid_name", ""), "worms_valid_AphiaID": w.get("valid_AphiaID", ""),
                "worms_status": w.get("status", ""), "worms_match_type": w.get("match_type", ""),
                "worms_class": w.get("class", ""), "worms_order": w.get("order", ""),
                "worms_family": w.get("family", "")}

    # --- assembly rows
    runs_by_asm = defaultdict(list)
    runs_by_tax = defaultdict(list)
    for r in runs:
        runs_by_tax[r["species_tax_id"]].append(r)
        for item in filter(None, r["linked_assemblies_derived"].split(";")):
            acc, ev = item.rsplit(":", 1)
            runs_by_asm[acc].append((r, ev))
    mito_by_tax = defaultdict(list)
    for m in mito_rows:
        for t in filter(None, m["species_tax_id"].split(";")):
            mito_by_tax[t].append(m)

    reps = {t: core.representative(v)["assembly_accession"] for t, v in by_tax.items()}
    asm_rows, qc = [], []
    rel = defaultdict(list)
    for a in asms:
        acc, t = a["assembly_accession"], a["tax_id"]
        linked = runs_by_asm.get(acc, [])
        row = {**taxo(t), **(worms_cols(tax.get(t, {}).get("species") or a["organism_name"]) if worms else {}),
               **a, "representative_derived": "Y" if reps[t] == acc else ""}
        found = set()
        for c in core.CATEGORIES:
            rs = [r for r, _ in linked if r["category_derived"] == c]
            if rs:
                found.add(c)
            row[f"{c}_runs"] = len(rs) or ""
            row[f"{c}_Gb"] = round(math.fsum(r["bases_Gb"] for r in rs), 2) or ""
            row[f"{c}_accessions"] = core.join(r["Run"] for r in rs)
        row["data_types_derived"] = "+".join(c for c in core.CATEGORIES if c in found)
        row["link_evidence_derived"] = core.join(sorted({ev for _, ev in linked}))
        sp_cats = {r["category_derived"] for r in runs_by_tax.get(t, [])}
        row["species_SRA_types_derived"] = "+".join(c for c in core.CATEGORIES if c in sp_cats)
        row["species_SRA_runs"] = len(runs_by_tax.get(t, [])) or ""
        mt_in = mito_asm.get(acc, [])
        linked_bs = {r.get("BioSample", "") for r, _ in linked}
        row.update(core.choose_mito(a["biosample"], linked_bs, mt_in, mito_by_tax.get(t, [])))
        for k, accs in row["_levels"].items():
            for m_acc in accs:
                rel[core._base(m_acc)].append(f"{k}:{acc}")
        asm_rows.append(row)

        # QC
        dec = core.declared_flags(a["sequencing_tech"])
        if "PacBio" in dec and not found & {"HiFi", "CLR", "PacBio_unknown"}:
            qc.append(("assembly", acc, a["organism_name"], f"sequencing_tech '{a['sequencing_tech']}' declares "
                       "PacBio but no linked PacBio run"))
        if "ONT" in dec and "ONT" not in found:
            qc.append(("assembly", acc, a["organism_name"], f"sequencing_tech '{a['sequencing_tech']}' declares "
                       "ONT but no linked ONT run"))
        if "Hi-C" in dec and "Hi-C" not in found:
            qc.append(("assembly", acc, a["organism_name"], f"sequencing_tech '{a['sequencing_tech']}' declares "
                       "Hi-C but no linked Hi-C run"))
        if not linked:
            qc.append(("assembly", acc, a["organism_name"], "no SRA run linked by BioSample/BioProject"
                       + (f" ({len(runs_by_tax.get(t, []))} species-level runs exist)" if runs_by_tax.get(t) else "")))
        if a["organelles_declared"].lower().find("mitochondrion") >= 0 and not mt_in:
            qc.append(("assembly", acc, a["organism_name"], "organelle_info declares a mitochondrion but the "
                       "sequence report has no mitochondrial molecule"))
        if row["mito_source_derived"].startswith("SUBSTITUTE"):
            qc.append(("mitogenome", acc, a["organism_name"], f"no mitogenome from the assembly's own sample; "
                       f"same-name substitute {row['mito_accession']} ({row['mito_note_derived']})"))
        if a["assembly_content_derived"] != "nuclear genome":
            qc.append(("assembly", acc, a["organism_name"], f"{a['assembly_content_derived']}: total length "
                       f"{a['Assembly Stats Total Sequence Length']} bp; not a nuclear genome"))
        if taxo(t)["ncbi_rank"] not in ("species", "subspecies"):
            qc.append(("taxonomy", t, a["organism_name"], f"NCBI rank is '{taxo(t)['ncbi_rank']}', not species"))
    if worms:
        for name, w in sorted(worms.items()):
            if not w:
                qc.append(("worms", "", name, "no WoRMS match"))
            elif w.get("status") != "accepted":
                qc.append(("worms", w.get("AphiaID", ""), name,
                           f"WoRMS status '{w.get('status')}' -> accepted name {w.get('valid_name')}"))
    pac_unknown = [r for r in runs if r["category_derived"] == "PacBio_unknown"]
    if pac_unknown:
        qc.append(("runs", "", "", f"{len(pac_unknown)} PacBio run(s) without HiFi/CLR evidence: "
                   + core.join((r["Run"] for r in pac_unknown), limit=20)))

    for m in mito_rows:
        m["relation_derived"] = core.join(rel.get(core._base(m["accession"]), [])) or "same species only"

    # --- species rows
    order = list(core.MITO_SOURCES.values())
    sp_rows = []
    for t, v in sorted(by_tax.items(), key=lambda kv: [taxo(kv[0])[k] for k in core.RANKS] + [kv[1][0]["organism_name"]]):
        rep = next(a for a in asm_rows if a["assembly_accession"] == reps[t])
        rs = runs_by_tax.get(t, [])
        linked_cats = {c for a in asm_rows if a["tax_id"] == t for c in a["data_types_derived"].split("+") if c}
        row = {**taxo(t), **(worms_cols(tax.get(t, {}).get("species") or v[0]["organism_name"]) if worms else {}),
               "organism_name": v[0]["organism_name"], "tax_id": t, "n_assemblies": len(v),
               "assemblies": core.join(a["assembly_accession"] for a in v),
               "representative_assembly": rep["assembly_accession"], "rep_paired_accession": rep["paired_accession"],
               "rep_assembly_level": rep["assembly_level"], "rep_total_length_Mb": rep["total_length_Mb"],
               "rep_contig_N50_kb": rep["contig_N50_kb"], "rep_scaffold_N50_kb": rep["scaffold_N50_kb"],
               "rep_sequencing_tech": rep["sequencing_tech"], "rep_release_date": rep["release_date"],
               "rep_submitter": rep["submitter"],
               "assembly_linked_types_derived": "+".join(c for c in core.CATEGORIES if c in linked_cats)}
        for c in core.CATEGORIES:
            cr = [r for r in rs if r["category_derived"] == c]
            row[f"species_{c}_runs"] = len(cr) or ""
            row[f"species_{c}_Gb"] = round(math.fsum(r["bases_Gb"] for r in cr), 2) or ""
        # species mitogenome: the best provenance among all its assemblies (representative first on ties)
        cand = sorted((a for a in asm_rows if a["tax_id"] == t and a["mito_accession"]),
                      key=lambda a: (order.index(a["mito_source_derived"]), a["assembly_accession"] != reps[t]))
        best = cand[0] if cand else {}
        row.update({k: best.get(k, "") for k in ("mito_accession", "mito_length_bp", "mito_source_derived",
                                                 "mito_note_derived")})
        row["mito_from_assembly"] = best.get("assembly_accession", "")
        row["species_mitogenomes_n"] = len(mito_by_tax.get(t, [])) or ""
        row["species_mitogenomes"] = core.join((m["accession"] for m in mito_by_tax.get(t, [])), limit=20)
        sp_rows.append(row)

    return asm_rows, sp_rows, runs, mito_rows, qc


def enrich(asm_rows, sp_rows, runs, sample_info, goat):
    """Same-individual evidence, category confidence, GoaT status; returns (goat_pipeline rows, qc rows)."""
    graph = core.sample_graph(sample_info.get("relations", {}))
    attrs = sample_info.get("attrs", {})
    asm_bs = {a["assembly_accession"]: a["biosample"] for a in asm_rows}
    per_asm = defaultdict(lambda: defaultdict(set))       # acc -> label -> categories
    for r in runs:
        r["category_confidence_derived"] = core.category_confidence(r["category_rule_derived"])
        ev = []
        for item in filter(None, r["linked_assemblies_derived"].split(";")):
            acc = item.rsplit(":", 1)[0]
            e = core.individual_evidence(r.get("BioSample", ""), asm_bs.get(acc, ""), graph, attrs)
            ev.append(f"{acc}: {e}")
            per_asm[acc][e.split(":")[0]].add(r["category_derived"])
        r["individual_derived"] = ";".join(ev)
    for a in asm_rows:
        d = per_asm.get(a["assembly_accession"], {})
        for key, col in (("same individual", "same_individual_types_derived"),
                         ("different individual", "other_individual_types_derived"),
                         ("unknown", "individual_unknown_types_derived")):
            a[col] = "+".join(c for c in core.CATEGORIES if c in d.get(key, set()))
    qc, pipeline = [], []
    if goat:
        status = {g["tax_id"]: g["goat_sequencing_status"] for g in goat["species"]}
        for srow in sp_rows:
            srow["goat_sequencing_status"] = status.get(srow["tax_id"], "")
        ours = {a["assembly_accession"] for a in asm_rows} | {a["paired_accession"] for a in asm_rows}
        for g in goat["assemblies"]:
            if g["accession"] and g["accession"] not in ours:
                qc.append(("crosscheck", g["accession"], g["scientific_name"], "in GoaT but not in our NCBI list"))
        have = {s["tax_id"] for s in sp_rows}
        pipeline = [g for g in goat["species"] if g["tax_id"] not in have and g["goat_sequencing_status"]
                    in ("sample_collected", "sample_acquired", "data_generation", "in_progress", "in_assembly")]
    return pipeline, qc


def write_outputs(out, taxon, taxid, asm_rows, sp_rows, runs, mito_rows, qc, args, counts, pipeline=()):
    date = f"{datetime.date.today():%Y%m%d}"
    head_tax = list(core.RANKS) + ["ncbi_rank"]
    worms_c = ["worms_valid_name", "worms_valid_AphiaID", "worms_status", "worms_match_type",
               "worms_class", "worms_order", "worms_family"] if not args.no_worms else []
    asm_cols = core.NCBI_COLS + head_tax + worms_c + [
        "representative_derived", "found_in_derived", "assembly_content_derived", "assembly_type", "refseq_category", "bioproject_parents", "tissue",
        "geo_loc_name", "assembly_method", "genome_coverage", "annotation_provider",
        "data_types_derived", "link_evidence_derived", "same_individual_types_derived",
        "other_individual_types_derived", "individual_unknown_types_derived"]
    asm_cols += [f"{c}_{s}" for c in core.CATEGORIES for s in ("runs", "Gb", "accessions")]
    asm_cols += ["species_SRA_types_derived", "species_SRA_runs", "organelles_declared", "mito_accession",
                 "mito_length_bp", "mito_source_derived", "mito_note_derived", "mito_same_sample_all",
                 "mito_substitute_n", "mito_substitute_candidates", "linked_assemblies",
                 "assembly_status"]
    # drop category columns that are empty everywhere (keeps the sheet readable)
    asm_cols = [c for c in asm_cols if not any(c.startswith(p + "_") for p in core.CATEGORIES)
                or any(r.get(c) not in ("", None) for r in asm_rows)]
    sp_cols = list(sp_rows[0].keys()) if sp_rows else []
    sp_cols = [c for c in sp_cols if not c.startswith("species_") or c in ("species_mitogenomes_n", "species_mitogenomes")
               or any(r.get(c) not in ("", None) for r in sp_rows)]
    run_cols = list(runs[0].keys()) if runs else ["Run"]
    mito_cols = ["species_tax_id", "accession", "organism", "tax_id", "length", "topology", "completeness",
                 "title", "biosample", "create_date", "source_db", "relation_derived"]
    qc_cols = ["scope", "id", "name", "message"]
    query = [("taxon", taxon), ("taxid", taxid), ("date", date),
             ("assemblies", "all versions" if args.all_versions else "current version only (GCA; paired GCF in paired_accession)"),
             ("sources", "NCBI Datasets v2 (primary; replaced the Assembly website in June 2024) + Entrez assembly db "
                         "(latest_genbank) as cross-check; union kept, see found_in_derived and QC 'crosscheck'"),
             ("worms", "off" if args.no_worms else "AphiaRecordsByMatchNames, first candidate, marine_only=false"),
             ("mitogenome records", f"species: nuccore mitochondrion[filter] AND {MITO_MIN}:{MITO_MAX}[SLEN]; "
                                    "plus every mitochondrion[filter] record of the assembly / raw-data BioSamples "
                                    "(any length)"),
             ("mitogenome rule", "one per assembly, in this order: " + " > ".join(core.MITO_SOURCES.values())
              + ". SUBSTITUTE = same species name only; NCBI identifications are not checked, verify before use."),
             ("link rule", "biosample = same BioSample; bioproject = same project or parent/child; "
                           "sibling_project = run project and assembly project share a direct parent "
                           "(e.g. Sanger/DToL species umbrella); species-level runs without link are only counted "
                           "in species_* columns"),
             ("category rule", "Hi-C: LibraryStrategy=Hi-C or Hi-C/Omni-C/Arima/Dovetail keyword; PacBio: "
                               "HiFi/CCS keyword -> HiFi, CLR/subread keyword -> CLR, Revio/Vega -> HiFi, "
                               "RS/RS II/Sequel -> CLR, else PacBio_unknown; ONT genomic/transcriptomic; "
                               "transcriptomic short -> RNA-seq; WGS short -> short_read_WGS; 10x/stLFR keyword -> "
                               "10x_linked. Keywords from experiment title, design, library name/protocol."),
             ("same individual", "run sample = assembly sample, or linked by ENA sample relation "
                                 "(derived_from/same_as/has_part, directly or via a common parent), or same ToLID; "
                                 "different ToLIDs = different individual; otherwise unknown"),
             ("category_confidence_derived", "high = explicit strategy/protocol/keyword or HiFi-only instrument; "
                                             "medium = pre-HiFi instrument; low = inferred (Sequel IIe) or unresolved"),
             ("GoaT", "goat.genomehubs.org: assembly cross-check + species in sequencing pipeline without an NCBI "
                      "assembly (sheet GoaT_pipeline)"),
             *[(k, v) for k, v in counts.items()]]
    prev = out / f"{taxon}_assemblies.tsv"
    if prev.exists():
        with open(prev, encoding="utf-8") as f:
            changes = core.change_report(list(csv.DictReader(f, delimiter="\t")), asm_rows)
        if not changes:
            changes = [{"change": "no change since the previous run"}]
    else:
        changes = [{"change": "first run (nothing to compare)"}]
    ch_cols = ["change", "assembly", "organism", "field", "old", "new"]
    pl_cols = ["class", "order", "family", "scientific_name", "tax_id", "goat_sequencing_status", "goat_projects"]
    xlsx = out / f"{taxon}_genome_inventory_{date}.xlsx"
    xlsx = write_xlsx(xlsx, {
        "Assemblies": (asm_cols, asm_rows), "Species": (sp_cols, sp_rows), "Runs": (run_cols, runs),
        "Mitogenomes": (mito_cols, mito_rows), "QC": (qc_cols, [dict(zip(qc_cols, q)) for q in qc]),
        "Changes": (ch_cols, changes), "GoaT_pipeline": (pl_cols, pipeline),
        "Query": (["item", "value"], [list(q) for q in query])})
    for name, cols, rows in (("assemblies", asm_cols, asm_rows), ("species", sp_cols, sp_rows),
                             ("runs", run_cols, runs), ("mitogenomes", mito_cols, mito_rows),
                             ("qc", qc_cols, [dict(zip(qc_cols, q)) for q in qc]),
                             ("changes", ch_cols, changes), ("goat_pipeline", pl_cols, pipeline)):
        with open(out / f"{taxon}_{name}.tsv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows({k: core.tsv_safe(v) for k, v in r.items()} for r in rows)
    return xlsx


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("taxon", help="taxon name (e.g. Annelida) or NCBI taxid")
    p.add_argument("-o", "--out", help="output directory (default: <taxon>_genome_inventory)")
    p.add_argument("--taxid", help="NCBI taxid when the name is ambiguous")
    p.add_argument("--all-versions", action="store_true", help="include superseded assembly versions")
    p.add_argument("--no-worms", action="store_true", help="skip WoRMS name matching")
    p.add_argument("--no-goat", action="store_true", help="skip the GoaT cross-check")
    p.add_argument("--update", action="store_true",
                   help="refetch the assembly lists; SRA/mitogenome/taxonomy data only for new species/samples")
    p.add_argument("--refresh", action="store_true", help="ignore all cached raw responses and refetch everything")
    p.add_argument("--email", help="NCBI e-mail (else $NCBI_EMAIL)")
    p.add_argument("--api-key", help="NCBI API key (else $NCBI_API_KEY)")
    args = p.parse_args(argv)

    taxon = args.taxon.strip().replace(" ", "_")
    out = Path(args.out or f"{taxon}_genome_inventory")
    cache = Cache(out / "raw", args.refresh)
    cl = Client(args.email, args.api_key)
    try:
        taxid = args.taxid or cache.get("taxid", lambda: resolve_taxid(cl, args.taxon))
        log(f"{args.taxon} -> taxid {taxid}")
        reps = cache.get("assemblies" + ("_all" if args.all_versions else ""),
                         lambda: fetch_assemblies(cl, taxid, args.all_versions), refresh=args.update)
        entrez = cache.get("entrez_assemblies" + ("_all" if args.all_versions else ""),
                           lambda: fetch_entrez_assemblies(cl, taxid, args.all_versions), refresh=args.update)
        ds_accs = {r["accession"] for r in reps} | {r.get("paired_accession", "") for r in reps}
        only_entrez = sorted({e["accession"] for e in entrez} - ds_accs)
        if only_entrez:
            log(f"{len(only_entrez)} assemblies only in Entrez -> fetching their Datasets reports by accession")
            added = cache.get("assemblies_added", lambda: fetch_assemblies_by_accession(cl, only_entrez),
                              refresh=args.update)
            have = {r["accession"] for r in reps}
            reps = reps + [r for r in added if r["accession"] not in have]
        if not reps:
            log("no genome assembly found")
            return 1
        asm_taxids = sorted({str(r["organism"]["tax_id"]) for r in reps})
        with_mt = sorted({r["accession"] for r in reps if any(
            "mitochondr" in (o.get("description") or "").lower() for o in r.get("organelle_info") or [])})
        log(f"{len(reps)} assemblies, {len(asm_taxids)} taxa, {len(with_mt)} declare a mitochondrion")
        mito_only = {r["accession"] for r in reps
                     if 0 < core.num(r.get("assembly_stats", {}).get("total_sequence_length")) < 1e6}
        ftp = {e["accession"]: e.get("ftp", "") for e in entrez}
        scan = fetch_mito_in_assembly(cl, reps, ftp, mito_only, out / "raw" / "mito_scan.json")
        mito_asm = {acc: [tuple(x) for x in v["mito"]] for acc, v in scan.items()}
        scan_err = {acc: v["method"] for acc, v in scan.items() if v["method"].startswith("ERROR")}
        log(f"  {sum(1 for v in mito_asm.values() if v)} assemblies contain a mitochondrion; "
            f"{len(scan_err)} could not be scanned")
        log("SRA runs for these taxa ...")
        sra = cache.get_incremental("sra", asm_taxids, lambda xs: fetch_sra(cl, xs),
                                    lambda a, b: {"runinfo": a["runinfo"] + b["runinfo"], "xml": {**a["xml"], **b["xml"]}})
        log(f"  {len(sra['runinfo'])} runs")
        run_taxids = {core.norm_taxid(r.get("TaxID")) for r in sra["runinfo"]} - {""}
        tax = cache.get_incremental("taxonomy", set(asm_taxids) | run_taxids, lambda xs: fetch_taxonomy(cl, xs),
                                    lambda a, b: {**a, **b})
        asm_bp = {r["assembly_info"].get("bioproject_accession", "") for r in reps}
        run_bp = {r.get("BioProject", "") for r in sra["runinfo"]}
        bp_parents = cache.get_incremental("bioprojects", (asm_bp | run_bp) - {""}, lambda xs: fetch_bioprojects(cl, xs),
                                           lambda a, b: {**a, **b})
        log("mitogenome records ...")
        mitos = cache.get_incremental("mitogenomes", asm_taxids, lambda xs: fetch_mitogenomes(cl, xs),
                                      lambda a, b: a + b)
        biosamples = ({(r["assembly_info"].get("biosample") or {}).get("accession", "") for r in reps}
                      | {r.get("BioSample", "") for r in sra["runinfo"]}) - {""}
        by_bs = cache.get_incremental("mito_by_biosample", biosamples, lambda xs: fetch_mito_by_biosample(cl, xs),
                                      lambda a, b: a + b)
        seen = {m["accession"] for m in mitos}
        mitos = mitos + [m for m in by_bs if m["accession"] not in seen]
        log(f"  {len(mitos)} mitochondrial records ({len(by_bs)} found by BioSample of {len(biosamples)})")
        worms = {}
        if not args.no_worms:
            names = {tax.get(t, {}).get("species") or "" for t in asm_taxids}
            names |= {str(r["organism"]["organism_name"]) for r in reps if not tax.get(str(r["organism"]["tax_id"]), {}).get("species")}
            worms = cache.get_incremental("worms", names - {""}, lambda xs: fetch_worms(cl, xs),
                                          lambda a, b: {**a, **b})
        log("ENA sample relations / ToLID ...")
        sample_info = cache.get_incremental(
            "sample_info", biosamples, lambda xs: fetch_sample_info(cl, xs),
            lambda a, b: {"relations": {**a["relations"], **b["relations"]}, "attrs": {**a["attrs"], **b["attrs"]}})
    except NcbiError as e:
        log(f"FAILED: {e}\n  rerun the same command to resume from {out / 'raw'}")
        return 1
    goat = {}
    if not args.no_goat:
        try:
            goat = cache.get("goat", lambda: fetch_goat(cl, taxid), refresh=args.update)
        except (NcbiError, ValueError, KeyError) as e:
            log(f"GoaT skipped: {e}"[:200])

    asm_rows, sp_rows, runs, mito_rows, qc = build(args.taxon, taxid, reps, mito_asm, tax, sra, bp_parents,
                                                   mitos, worms, args, entrez, set(only_entrez))
    qc += core.crosscheck(asm_rows, entrez)
    pipeline, qc_goat = enrich(asm_rows, sp_rows, runs, sample_info, goat)
    qc += qc_goat
    qc += [("assembly", acc, "", f"sequence list not scanned for a mitochondrion ({msg}); rerun to retry")
           for acc, msg in sorted(scan_err.items())]
    counts = {
        "assemblies (rows)": len(asm_rows), "species/taxa": len(sp_rows),
        "SRA runs of these taxa": len({r["Run"] for r in runs}),
        "runs linked to an assembly": len({r["Run"] for r in runs if r["linked_assemblies_derived"]}),
        "mitochondrial records": len(mito_rows), "QC items": len(qc),
        "GoaT species in pipeline without NCBI assembly": len(pipeline)}
    xlsx = write_outputs(out, taxon, taxid, asm_rows, sp_rows, runs, mito_rows, qc, args, counts, pipeline)
    types = defaultdict(int)
    for a in asm_rows:
        for c in filter(None, a["data_types_derived"].split("+")):
            types[c] += 1
    log("DONE")
    for k, v in counts.items():
        print(f"  {k}: {v}")
    print("  assemblies with linked data type: " + ", ".join(f"{c} {types[c]}" for c in core.CATEGORIES if types[c]))
    src = defaultdict(int)
    for a in asm_rows:
        src[a["mito_source_derived"] or "none"] += 1
    print("  mitogenome per assembly: " + "; ".join(f"{k} {v}" for k, v in src.items()))
    ind = defaultdict(int)
    for r in runs:
        for item in filter(None, r["individual_derived"].split(";")):
            ind[item.split(": ", 1)[1].split(":")[0]] += 1
    print("  linked run-assembly pairs by individual: " + "; ".join(f"{k} {v}" for k, v in ind.items()))
    print(f"  -> {xlsx}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
