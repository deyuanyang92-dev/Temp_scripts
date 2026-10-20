"""Build the Metazoa-wide NC_ (RefSeq) mitogenome benchmark set.

Per phylum: NCBI nuccore "<phylum>[Organism] AND mitochondrion[Filter] AND complete genome[Title] AND
refseq[Filter]", NC_ only. A random candidate pool (seeded) is fetched, ranked by NCBI Taxonomy, and
`pick` takes <= 2 records per order, tier A (transcriptome / long-read evidence in the source INSDC
record) first, up to the phylum quota. dev/test are split by genus (no genus in both).

Writes <out>/gb/<acc>.gb (NCBI record), <out>/bare/<acc>/<acc>.gbf (source feature only) and the
record list (default benchmark/sets/mito_bench_<YYYY-MM>.tsv). Re-running reuses downloaded files.

Usage: python3 mito_sample.py --out projects/mito_benchmark [--seed 20261020] [--tsv PATH] [--phyla A,B]
"""
import argparse
import csv
import datetime
import io
import random
import re
import sys
import time
import urllib.parse
import urllib.request
import defusedxml.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

BIG = 12
SMALL = 6
PHYLA: Dict[str, int] = {p: BIG for p in ("Chordata", "Arthropoda", "Mollusca", "Annelida", "Cnidaria",
                                           "Echinodermata", "Nematoda", "Platyhelminthes")}
PHYLA.update({p: SMALL for p in ("Porifera", "Placozoa", "Ctenophora", "Xenacoelomorpha", "Hemichordata",
                                 "Tardigrada", "Onychophora", "Nematomorpha", "Priapulida", "Kinorhyncha",
                                 "Loricifera", "Nemertea", "Brachiopoda", "Phoronida", "Bryozoa", "Entoprocta",
                                 "Cycliophora", "Rotifera", "Acanthocephala", "Gastrotricha", "Gnathostomulida",
                                 "Chaetognatha", "Orthonectida", "Dicyemida")})
# MitoZ 3.6 `annotate --clade` (batch_mitoz.py MITOZ_CLADES); None = not supported by MitoZ
MITOZ_CLADE: Dict[str, Optional[str]] = {p: None for p in PHYLA}
MITOZ_CLADE.update({"Chordata": "Chordata", "Arthropoda": "Arthropoda", "Echinodermata": "Echinodermata",
                    "Annelida": "Annelida-segmented-worms", "Bryozoa": "Bryozoa", "Mollusca": "Mollusca",
                    "Nematoda": "Nematoda", "Nemertea": "Nemertea-ribbon-worms", "Porifera": "Porifera-sponges"})
# fallback genetic code when a record has no CDS transl_table (flagged code_assumed)
DEFAULT_CODE = defaultdict(lambda: 5, {"Chordata": 2, "Echinodermata": 9, "Hemichordata": 9, "Platyhelminthes": 9,
                                       "Cnidaria": 4, "Porifera": 4, "Placozoa": 4, "Ctenophora": 4})
COLUMNS = ["accession", "source_accession", "phylum", "class", "order", "family", "genus", "species", "genetic_code",
           "code_assumed", "topology", "length", "evidence_tier", "annotation_method", "refseq_created",
           "structure_flags", "mitoz_clade", "split", "retrieved"]
ORDER_CAP = 2
POOL = 8                       # candidates fetched per quota slot
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOLS = ("MITOS", "MitoZ", "MFannot", "DOGMA", "GeSeq", "MitoAnnotator", "MitoFinder", "MitoMaker")


# ── record properties ────────────────────────────────────────────────────────
def majority_code(rec: SeqRecord, phylum: str) -> Tuple[int, bool]:
    """Most frequent CDS transl_table, or (phylum default, True) when no CDS gives one."""
    codes = [f.qualifiers["transl_table"][0] for f in rec.features if f.type == "CDS" and "transl_table" in f.qualifiers]
    if not codes:
        return DEFAULT_CODE[phylum], True
    return int(Counter(codes).most_common(1)[0][0]), False


def _text(rec: Optional[SeqRecord]) -> str:
    if rec is None:
        return ""
    ann: Dict[str, Any] = rec.annotations
    sc: Dict[str, Dict[str, str]] = ann.get("structured_comment") or {}
    parts = [str(ann.get("comment", ""))]
    parts += [f"{k} :: {v}" for block in sc.values() for k, v in block.items()]
    parts += [getattr(r, "title", "") for r in ann.get("references", [])]
    return "\n".join(parts)


def evidence_tier(source_rec: Optional[SeqRecord]) -> str:
    """A-longread (PacBio / Nanopore), A-transcriptome (RNA-seq, transcriptome, Iso-Seq, cDNA) or B."""
    t = _text(source_rec)
    if re.search(r"pacbio|nanopore|oxford|\bont\b|minion|promethion|hifi|sequel", t, re.I):
        return "A-longread"
    if re.search(r"rna-?seq|transcriptom|iso-?seq|\bcdna\b", t, re.I):
        return "A-transcriptome"
    return "B"


def annotation_method(rec: SeqRecord, source_rec: Optional[SeqRecord]) -> str:
    """Annotation tool named in the RefSeq / source record, 'other:<text>' or 'unknown'."""
    t = _text(rec) + "\n" + _text(source_rec)
    for tool in TOOLS:
        if re.search(rf"\b{tool}", t, re.I):
            return tool
    m = re.search(r"Annotation (?:Method|Pipeline) :: ([^\n]+)", t)
    return f"other:{m.group(1).strip()}" if m else "unknown"


def source_accession(rec: SeqRecord) -> str:
    m = re.search(r"(?:identical to|derived from)\s+([A-Z]{1,2}\d{5,8}(?:\.\d+)?)", str(rec.annotations.get("comment", "")))
    return m.group(1).rstrip(".") if m else ""


def structure_flags(rec: SeqRecord) -> List[str]:
    out = []
    if rec.annotations.get("topology") == "linear":
        out.append("linear")
    names = {q.lower() for f in rec.features for q in f.qualifiers.get("gene", [])}
    if not names & {"cox1", "coi", "co1", "coxi"}:
        out.append("no_cox1")
    if re.search(r"chromosome|minicircle|segment|fragment", str(rec.description), re.I):
        out.append("multi_chromosome")
    return out


# ── selection ────────────────────────────────────────────────────────────────
def pick(rows: List[dict], quota: int, seed: int) -> List[dict]:
    """<= ORDER_CAP per order (tier A first, then random); fill to quota from the rest (tier A first)."""
    rng = random.Random(seed)
    pool = rows[:]
    rng.shuffle(pool)
    pool.sort(key=lambda r: not r["evidence_tier"].startswith("A"))      # stable: A first, random within
    got: List[dict] = []
    rest: List[dict] = []
    per: Counter = Counter()
    for r in pool:
        if len(got) < quota and per[r["order"]] < ORDER_CAP:
            got.append(r)
            per[r["order"]] += 1
        else:
            rest.append(r)
    return got + rest[:max(0, quota - len(got))]


def split_by_genus(rows: List[dict], seed: int) -> None:
    """Set row['split'] to dev/test; genera alternate within each phylum after a seeded shuffle."""
    rng = random.Random(seed)
    by_phylum = defaultdict(list)
    for r in rows:
        by_phylum[r["phylum"]].append(r)
    for phylum in sorted(by_phylum):
        genera = sorted({r["genus"] for r in by_phylum[phylum]})
        rng.shuffle(genera)
        side = {g: ("dev", "test")[i % 2] for i, g in enumerate(genera)}
        for r in by_phylum[phylum]:
            r["split"] = side[r["genus"]]


def ranks_from_lineage_xml(xml: str) -> Dict[str, Dict[str, str]]:
    """taxid -> {phylum, class, order, family, genus} from NCBI Taxonomy efetch XML."""
    out = {}
    for t in ET.fromstring(xml).findall("Taxon"):
        ranks = {x.findtext("Rank"): x.findtext("ScientificName") for x in t.findall("LineageEx/Taxon")}
        out[t.findtext("TaxId") or ""] = {k: ranks.get(k) or "" for k in ("phylum", "class", "order", "family", "genus")}
    return out


# ── NCBI ─────────────────────────────────────────────────────────────────────
def eutil(name: str, **params) -> str:
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(5):
        try:
            time.sleep(0.34)
            with urllib.request.urlopen(f"{EUTILS}/{name}.fcgi", data=data, timeout=120) as r:
                return r.read().decode()
        except Exception as e:                       # network: retry, then give up loudly
            if attempt == 4:
                raise RuntimeError(f"{name} failed: {e}") from e
            time.sleep(3 * (attempt + 1))
    return ""


def search_ids(phylum: str) -> List[str]:
    term = f"{phylum}[Organism] AND mitochondrion[Filter] AND complete genome[Title] AND refseq[Filter]"
    xml = eutil("esearch", db="nuccore", term=term, retmax=100000)
    return [x.text or "" for x in ET.fromstring(xml).findall("IdList/Id")]


def summaries(ids: List[str]) -> List[dict]:
    out = []
    for i in range(0, len(ids), 400):
        root = ET.fromstring(eutil("esummary", db="nuccore", id=",".join(ids[i:i + 400]), version="2.0"))
        for d in root.iter("DocumentSummary"):
            out.append({"accession": d.findtext("AccessionVersion") or "", "taxid": d.findtext("TaxId") or "",
                        "created": (d.findtext("CreateDate") or "").replace("/", "-"), "title": d.findtext("Title") or ""})
    return out


def fetch_gb(accs: List[str], gb_dir: Path) -> Dict[str, SeqRecord]:
    need = [a for a in accs if not (gb_dir / f"{a}.gb").exists()]
    for i in range(0, len(need), 100):
        text = eutil("efetch", db="nuccore", id=",".join(need[i:i + 100]), rettype="gbwithparts", retmode="text")
        for r in SeqIO.parse(io.StringIO(text), "genbank"):
            SeqIO.write(r, str(gb_dir / f"{r.id}.gb"), "genbank")
    return {a: SeqIO.read(str(gb_dir / f"{a}.gb"), "genbank") for a in accs if (gb_dir / f"{a}.gb").exists()}


def taxonomy(taxids: List[str]) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    for i in range(0, len(taxids), 200):
        out.update(ranks_from_lineage_xml(eutil("efetch", db="taxonomy", id=",".join(taxids[i:i + 200]), retmode="xml")))
    return out


def write_bare(rec: SeqRecord, bare_dir: Path) -> None:
    d = bare_dir / str(rec.id)
    d.mkdir(parents=True, exist_ok=True)
    b = rec[:]
    b.annotations = dict(rec.annotations)
    b.features = [f for f in rec.features if f.type == "source"]
    SeqIO.write(b, str(d / f"{rec.id}.gbf"), "genbank")


# ── driver ───────────────────────────────────────────────────────────────────
def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True, help="data root (gb/, bare/, src/ are created inside)")
    p.add_argument("--seed", type=int, default=20261020)
    p.add_argument("--tsv", help="record list (default benchmark/sets/mito_bench_<YYYY-MM>.tsv)")
    p.add_argument("--phyla", help="comma-separated subset of phyla (default: all)")
    a = p.parse_args(argv)
    today = datetime.date.today().isoformat()
    out = Path(a.out)
    gb_dir, src_dir, bare_dir = out / "gb", out / "src", out / "bare"
    for d in (gb_dir, src_dir, bare_dir):
        d.mkdir(parents=True, exist_ok=True)
    tsv = Path(a.tsv) if a.tsv else Path(__file__).resolve().parent / "sets" / f"mito_bench_{today[:7]}.tsv"
    tsv.parent.mkdir(parents=True, exist_ok=True)
    phyla = a.phyla.split(",") if a.phyla else list(PHYLA)

    chosen = []
    for phylum in phyla:
        summ = [s for s in summaries(search_ids(phylum)) if s["accession"].startswith("NC_")]
        rng = random.Random(f"{a.seed}-{phylum}")
        rng.shuffle(summ)
        cand = sorted(summ[:PHYLA[phylum] * POOL], key=lambda s: s["accession"])
        if not cand:
            print(f"[sample] {phylum}: no NC_ records", file=sys.stderr)
            continue
        tax = taxonomy(sorted({s["taxid"] for s in cand}))
        recs = fetch_gb([s["accession"] for s in cand], gb_dir)
        srcs = fetch_gb(sorted({source_accession(r) for r in recs.values()} - {""}), src_dir)
        rows = []
        for s in cand:
            r = recs.get(s["accession"])
            if r is None:
                continue
            src_acc = source_accession(r)
            src: Optional[SeqRecord] = srcs.get(src_acc) or next((v for k, v in srcs.items() if k.split(".")[0] == src_acc.split(".")[0]), None)
            rk: Dict[str, str] = tax.get(s["taxid"], {})
            code, assumed = majority_code(r, phylum)
            species = str(r.annotations.get("organism", ""))
            rows.append({"accession": r.id, "source_accession": src_acc, "phylum": phylum,
                         "class": rk.get("class", ""), "order": rk.get("order", "") or f"?{rk.get('class', '')}",
                         "family": rk.get("family", ""), "genus": rk.get("genus", "") or species.split(" ")[0],
                         "species": species, "genetic_code": code, "code_assumed": "yes" if assumed else "",
                         "topology": r.annotations.get("topology", ""), "length": len(r),
                         "evidence_tier": evidence_tier(src) if src is not None else evidence_tier(r),
                         "annotation_method": annotation_method(r, src), "refseq_created": s["created"],
                         "structure_flags": ",".join(structure_flags(r)), "mitoz_clade": MITOZ_CLADE[phylum] or "",
                         "split": "", "retrieved": today})
        got = pick(rows, PHYLA[phylum], a.seed)
        print(f"[sample] {phylum}: {len(summ)} NC_, {len(rows)} candidates, {len(got)} chosen "
              f"(tier A {sum(r['evidence_tier'] != 'B' for r in got)})", file=sys.stderr)
        chosen += got
    split_by_genus(chosen, a.seed)
    for row in chosen:
        write_bare(SeqIO.read(str(gb_dir / f"{row['accession']}.gb"), "genbank"), bare_dir)
    with tsv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(sorted(chosen, key=lambda r: (r["phylum"], r["accession"])))
    print(f"[sample] {len(chosen)} records -> {tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
