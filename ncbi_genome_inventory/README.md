# ncbi_genome_inventory

One command → an Excel inventory of every NCBI whole-genome assembly under a taxon, with the raw data behind
each assembly (HiFi, Hi-C, RNA-seq, ONT …), the mitogenomes, and the taxonomy.

一条命令整理某类群（如环节动物）在 NCBI 的全部全基因组：登录号、所用原始数据（HiFi / Hi-C / RNA / ONT …）、
产出的线粒体、分类阶元。替代以前"网页勾选导出 + kingfisher + WoRMS 网页匹配 + 手工合并"的流程。

## Install / run

Python ≥ 3.8, **standard library only** (no pip install, no `datasets`/`edirect`/`kingfisher`).

```bash
python3 ncbi_genome_inventory/genome_inventory.py Annelida -o Annelida_genome_inventory
python3 ncbi_genome_inventory/genome_inventory.py Polynoidae --no-worms     # skip WoRMS
python3 ncbi_genome_inventory/genome_inventory.py 6340                      # taxid also works
python3 ncbi_genome_inventory/genome_inventory.py Annelida --update         # new NCBI assemblies (fast)
python3 ncbi_genome_inventory/genome_inventory.py Annelida --refresh        # refetch everything
```

- Raw responses are cached in `<out>/raw/`. If a step fails (network), rerun the same command: finished steps
  are read from the cache. `--update` refetches the assembly lists and fetches SRA / mitogenome / taxonomy
  data only for species, BioSamples and BioProjects not seen before (caches are keyed by these items);
  `--refresh` refetches everything. Annelida: full run ≈ 11 min, cached rerun ≈ 2 s.
- Optional: `export NCBI_API_KEY=...` (and `NCBI_EMAIL`) — about 3× faster E-utilities.
- `--all-versions` also lists superseded assembly versions (default: current versions only).
- Time (no key): Polynoidae 27 assemblies ≈ 1.5 min.

## Web interface

```bash
python3 ncbi_genome_inventory/webapp.py --data-dir /mnt/p/genome_inventory     # open http://127.0.0.1:8765
```

Type a taxon, press **Run**: the page runs the same `genome_inventory.py` (CLI and web always agree), shows the
log, then the tables (Assemblies / Species / Mitogenomes / QC / Runs) with filter, sort, "key columns only",
substitute mitogenomes highlighted, and download links for the Excel and TSV files. Existing result folders in
`--data-dir` are listed and can be opened. Listens on 127.0.0.1 only.

## Data sources (several, cross-checked)

| What | Method 1 | Method 2 / cross-check |
|---|---|---|
| assembly list + metadata | NCBI Datasets v2 REST (NCBI retired the Assembly/Genome websites in June 2024 and moved them to Datasets) | Entrez `assembly` db (`latest_genbank`; E-utilities keep working). Union kept; assemblies missing from Datasets' taxon query are fetched by accession; differences → QC `crosscheck`; column `found_in_derived` |
| sequences inside each assembly (mitochondrion) | FTP `*_assembly_report.txt` (one request, all sequence names) | Datasets `sequence_reports` if FTP fails; a failure is listed in QC and retried on the next run |
| raw data | SRA RunInfo + SRA XML (E-utilities) | ENA sample relations + ToLID for same-individual evidence |
| assemblies (third source) | — | GoaT assemblies and species sequencing status |
| mitogenome records | nuccore by species (10–100 kb) | nuccore by BioSample (any length) |
| taxonomy | NCBI Taxonomy | WoRMS (accepted name, class/order/family) |

## Output

`<out>/<Taxon>_genome_inventory_<YYYYMMDD>.xlsx` (+ the same tables as TSV):

| Sheet | One row per | Content |
|---|---|---|
| Assemblies | assembly (GCA; paired GCF in `Assembly Paired Assembly Accession`) | **first 41 columns = the NCBI Datasets web table / `dataformat tsv genome` columns, same names and order as the manual export** (Assembly Accession … CheckM contamination); then NCBI phylum–species, WoRMS accepted name/class/order/family, `representative_derived`, linked raw data per type (`HiFi_runs`, `HiFi_Gb`, `HiFi_accessions`, `Hi-C_…`, `RNA-seq_…` …), `data_types_derived` (e.g. `HiFi+Hi-C+RNA-seq`), `link_evidence_derived`, `mito_accession` + `mito_source_derived` (same sample first, see below) |
| Species | taxid | representative assembly (RefSeq reference > level > contig N50 > newest; alternate haplotypes never chosen), all assemblies, data types linked to assemblies, all SRA data of the species by type, mitogenomes |
| Runs | SRA run of those species | RunInfo fields, experiment title / library protocol, `category_derived` + `category_rule_derived`, `linked_assemblies_derived` (`GCA_…:evidence`) |
| Mitogenomes | nuccore mitochondrial record (species: 10–100 kb; assembly/raw-data BioSamples: any length) | `relation_derived`: `in_assembly:GCA_…` / `same_biosample:GCA_…` / `linked_run_biosample:GCA_…` / `substitute:GCA_…` / `same species only` |
| QC | issue | declared technology without linked run, assemblies with no public raw data, mitochondrion declared but missing, non-species taxa, WoRMS unmatched / unaccepted names, PacBio runs of unknown type |
| Query | — | parameters, rules, counts |

## How raw data are linked to an assembly

Only runs of the same species (taxid, or a subspecies below it) are considered, then:

1. `biosample` — run BioSample = assembly BioSample (same specimen);
2. `bioproject` — same BioProject, or one is the parent of the other;
3. `sibling_project` — run project and assembly project share a direct parent. This is how Sanger/DToL,
   ERGA etc. submit: e.g. *Acholoe squamosa* GCA_949317995.1 is in PRJEB60118, its HiFi and Hi-C reads are in
   PRJEB59385 from **other BioSamples**; both sit under the species umbrella PRJEB59386. Matching by BioSample
   alone would miss all of them (tested in `tests/test_core.py`).

Runs of the species with no such link are not attached to an assembly; they are counted in `species_*` columns.

## Mitogenome: same sample first (`mito_source_derived`)

One mitogenome per assembly, taken in this order — a same-name record is used **only** when no same-sample
record exists, and is flagged, because NCBI identifications are not reliable:

| `mito_source_derived` | Meaning |
|---|---|
| `same sample: in assembly` | mitochondrial molecule of the assembly itself (also detected by sequence name, since NCBI sometimes labels it `Chromosome`; a mitogenome-only assembly counts as its own mitogenome) |
| `same sample: same BioSample` | separate nuccore record with the assembly's BioSample (any length; < 10 kb noted as partial) |
| `same project sample: BioSample of linked raw data` | record from the BioSample of the HiFi/Hi-C/… runs linked to the assembly |
| `SUBSTITUTE: same species name, different sample - verify identity` | only the species name matches; listed in QC. Among candidates: INSDC original over RefSeq `NC_` copy > complete > circular > newest. All candidates are in `mito_substitute_candidates` |

Annelida 2026-10-28: 106 assemblies in-assembly, 65 same BioSample (mostly alternate haplotypes sharing the
primary's sample), 44 substitutes, 36 none.

## Data type rules (`category_derived`)

First match wins; evidence = LibraryStrategy/Source/Platform/Model + keywords in experiment title, design
description, library name and library construction protocol (never the study abstract):

| Category | Rule |
|---|---|
| Hi-C | LibraryStrategy = Hi-C, or keyword Hi-C / Omni-C / Arima / Dovetail / Phase Genomics / proximity ligation (not for RNA) |
| HiFi | PacBio genomic + keyword HiFi / CCS; or Revio / Vega; or Sequel IIe without keyword (*inferred*, rule says so) |
| CLR | PacBio genomic + keyword CLR / subread; or RS / RS II / Sequel (I) |
| PacBio_unknown | PacBio Sequel II without keyword (listed in QC) |
| ONT / ONT_RNA | Oxford Nanopore genomic / transcriptomic |
| Iso-Seq | PacBio transcriptomic |
| RNA-seq | short-read transcriptomic |
| short_read_WGS | WGS on Illumina / BGI / DNBSEQ … |
| 10x_linked | keyword 10x / Chromium / linked-read / stLFR / TELL-seq |
| Other | metagenomic, AMPLICON, RAD-Seq, … (`category_rule_derived` shows the strategy) |

Annelida check (317 PacBio experiments): almost all carry an explicit `PacBio - HiFi` protocol or HiFi/CCS/
subread keyword.

## Same individual? (`individual_derived`)

Project links say *where* reads were submitted, not *whose* reads they are. For every run linked to an assembly
the run's BioSample is compared with the assembly's BioSample:

- `same individual: same BioSample`;
- `same individual: BioSample relation` — ENA `related_sample_accession` (`derived_from` / `same_as` /
  `has_part`), directly or through a common parent sample. Live example: *Acholoe squamosa* HiFi sample
  SAMEA13854475 and Hi-C sample SAMEA13854477 are both `derived_from` the assembly sample SAMEA13854387;
- `same individual: same ToLID` / `different individual: ToLID x vs y` (a ToLID names one individual);
- `unknown`.

Assemblies get `same_individual_types_derived` / `other_individual_types_derived` /
`individual_unknown_types_derived`, e.g. `HiFi+Hi-C` from the reference individual, `RNA-seq` from another.

## Other sheets

- **Changes** — compared with the previous run in the same folder: new assemblies, assemblies no longer listed,
  changed level / N50 / BioSample / data types / mitogenome (idea from genome_updater's update report).
- **GoaT_pipeline** — species that GoaT lists as sample collected / acquired / data generation / in progress /
  in assembly but that have no NCBI assembly yet; GoaT assemblies missing from the NCBI list go to QC.
- `category_confidence_derived` on runs: high / medium / low (see the data type rules).

## Comparison with existing tools (2026-10-28)

No open-source tool does the whole chain *taxon → every assembly → its raw data (by sample/project) → same
individual → mitogenome from the same sample*. Closest:

| Tool | What it covers | Gap vs this tool |
|---|---|---|
| [GoaT](https://goat.genomehubs.org) / goat-cli | per-taxon assemblies, sequencing status, `mitochondrion_accession` | no per-assembly raw-data linkage or data types; used here as a cross-check + pipeline list |
| NCBI `datasets` / ncbi-genome-download / genome_updater | assembly metadata and downloads, update reports | no SRA linkage, no mitogenome provenance |
| kingfisher / pysradb / ffq / fastq-dl | SRA/ENA metadata and downloads for given accessions or projects | not taxon → assembly oriented; no HiFi/CLR/Hi-C classification against assemblies |
| Darwin Tree of Life portal | DToL species: runs, assemblies, ToLIDs, genome notes | DToL only, per species |

## Notes

- Values from NCBI are never changed; derived columns end with `_derived`.
- Some new assemblies (e.g. 2026 Sanger releases) have no public raw reads yet in NCBI SRA or ENA — QC lists
  them; that is not a linking failure.
- Tests (offline): `python3 -m pytest -q ncbi_genome_inventory/tests`.
- `ncbi.py` / `xlsx.py` are adapted from `yzz_taxonomy_workflow/modules/sra` (same author); the yzz `sra`
  module is the tool for a full SRA inventory of a taxon (all runs, shallow vs deep WGS).
