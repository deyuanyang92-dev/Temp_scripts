#!/usr/bin/env python3
"""fetch_rrna_refs.py - reference mitogenomes per phylum for the rRNA HMM library (mito_rrna_hmm.py build).

Reusable form of knowledge/experiments/rrna_ref/fetch_refs.py (which has the benchmark set and the dev
phyla hard-coded). Per phylum, NCBI nuccore is queried with
    <phylum>[Organism] AND mitochondrion[Filter] AND (complete genome[Title] OR complete sequence[Title])
    AND 5000:60000[SLEN]
a random sample (seeded, per phylum) of --pool of the first --retmax hits is downloaded (gbwithparts), records
of the exclusion set are removed and the first --per-phylum remaining ones are kept.

Output in -o:
  <phylum>.gb          GenBank records; the file stem is the phylum name that `mito_rrna_hmm.py build` needs
  refs_<phylum>.txt    accession.version of the kept records, one per line (audit trail of the library)
  fetch_log.tsv        phylum, query, date, hits, fetched, excluded, kept, seed (appended on every run)
An existing non-empty <phylum>.gb is not downloaded again (its accession list is rewritten).

--set TSV   records to leave out: every accession in the columns `accession` and `source_accession`
            (version ignored), e.g. benchmark/sets/mito_bench_2026-10.tsv; with --exclude-genera also every
            record whose organism genus is listed in the column `genus` (the stricter library of the
            2026-10-21 A/B experiment). Optional.

Usage:
  python3 scripts/fetch_rrna_refs.py -o refs --phyla Ctenophora Placozoa --set benchmark/sets/mito_bench_2026-10.tsv
  python3 scripts/fetch_rrna_refs.py -o refs --all-metazoa --per-phylum 150 --seed 7
  python3 mito_rrna_hmm.py build --gb refs/*.gb -o models/rrna_hmm
Needs network access (NCBI E-utilities, ~0.5 s between requests; large phyla take minutes).
Limit: esearch returns at most --retmax ids, so for phyla with more hits the sample is random only within
the first --retmax hits of the NCBI ordering.
"""
import argparse
import csv
import datetime
import json
import random
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import AbstractSet, Callable, Dict, List, Optional, Set, Tuple

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
QUERY = ("{}[Organism] AND mitochondrion[Filter] AND (complete genome[Title] OR complete sequence[Title]) "
         "AND 5000:60000[SLEN]")
# Metazoan phyla of the pipeline: copy of batch_mitoz.METAZOAN_PHYLA (kept separate on purpose; a test compares them)
METAZOAN_PHYLA = (
    "Chordata", "Arthropoda", "Echinodermata", "Annelida", "Bryozoa", "Mollusca", "Nematoda", "Nemertea", "Porifera",
    "Cnidaria", "Platyhelminthes", "Chaetognatha", "Hemichordata", "Xenacoelomorpha", "Ctenophora", "Placozoa",
    "Tardigrada", "Onychophora", "Nematomorpha", "Priapulida", "Kinorhyncha", "Loricifera", "Brachiopoda", "Phoronida",
    "Entoprocta", "Cycliophora", "Rotifera", "Acanthocephala", "Gastrotricha", "Gnathostomulida", "Orthonectida",
    "Dicyemida", "Micrognathozoa")
Get = Callable[[str], str]


def ncbi_get(url: str) -> str:
    for i in range(5):
        try:
            return urllib.request.urlopen(url, timeout=90).read().decode()
        except Exception:
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"NCBI request failed after 5 tries: {url}")


def load_set(path: Path) -> Tuple[Set[str], Set[str]]:
    """(accessions without version, genera) to leave out, from a TSV with accession / source_accession / genus."""
    acc: Set[str] = set()
    gen: Set[str] = set()
    with open(path, newline="") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            acc |= {r[c].strip().split(".")[0] for c in ("accession", "source_accession") if (r.get(c) or "").strip()}
            if (r.get("genus") or "").strip():
                gen.add(r["genus"].strip())
    return acc, gen


def accession_of(record: str) -> str:
    """accession.version of a GenBank record text (VERSION line)."""
    m = re.search(r"^VERSION\s+(\S+)", record, re.M)
    return m.group(1) if m else ""


def genus_of(record: str) -> str:
    m = re.search(r"^\s+ORGANISM\s+(\S+)", record, re.M)
    return m.group(1) if m else ""


def fetch_phylum(phylum: str, get: Get, per: int = 150, seed: int = 7, pool: int = 400, retmax: int = 5000,
                 exclude_acc: AbstractSet[str] = frozenset(), exclude_gen: AbstractSet[str] = frozenset(),
                 pause: float = 0.5) -> Tuple[List[str], dict]:
    """(GenBank record texts kept, log row) for one phylum; get(url) -> text is the E-utilities transport."""
    query = QUERY.format(phylum)
    js = json.loads(get(EUTILS + "esearch.fcgi?db=nuccore&retmax=%d&retmode=json&term=%s" % (retmax, urllib.parse.quote(query))))
    ids = list(js["esearchresult"]["idlist"])
    time.sleep(pause)
    random.Random(f"{seed}:{phylum}").shuffle(ids)
    ids = ids[:pool]
    recs: List[str] = []
    for i in range(0, len(ids), 100):
        txt = get(EUTILS + "efetch.fcgi?db=nuccore&rettype=gbwithparts&retmode=text&id=" + ",".join(ids[i:i + 100]))
        recs += [r + "//\n" for r in txt.split("//\n") if r.strip()]
        time.sleep(pause)
    keep = [r for r in recs if accession_of(r).split(".")[0] not in exclude_acc and genus_of(r) not in exclude_gen]
    return keep[:per], {"phylum": phylum, "query": query, "date": datetime.date.today().isoformat(),
                        "n_hits": int(js["esearchresult"]["count"]), "n_fetched": len(recs),
                        "n_excluded": len(recs) - len(keep), "n_kept": min(per, len(keep)), "seed": seed}


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", type=Path, required=True, help="output directory (<phylum>.gb, refs_<phylum>.txt)")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--phyla", nargs="+", help="phylum names (NCBI organism)")
    g.add_argument("--all-metazoa", action="store_true", help="all metazoan phyla of the pipeline")
    ap.add_argument("--set", type=Path, help="TSV with accession / source_accession (/ genus) to leave out")
    ap.add_argument("--exclude-genera", action="store_true", help="also leave out the genera of --set")
    ap.add_argument("--per-phylum", type=int, default=150, help="records kept per phylum (default 150)")
    ap.add_argument("--pool", type=int, default=400, help="random hits downloaded per phylum before exclusion")
    ap.add_argument("--retmax", type=int, default=5000, help="ids taken from the NCBI search (default 5000)")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--pause", type=float, default=0.5, help="seconds between NCBI requests")
    a = ap.parse_args(argv)
    acc, gen = load_set(a.set) if a.set else (set(), set())
    if not a.exclude_genera:
        gen = set()
    a.out.mkdir(parents=True, exist_ok=True)
    log: List[Dict[str, object]] = []
    for p in (list(METAZOAN_PHYLA) if a.all_metazoa else a.phyla):
        gb = a.out / f"{p}.gb"
        if gb.exists() and gb.stat().st_size:
            texts = [r + "//\n" for r in gb.read_text().split("//\n") if r.strip()]
            print(f"{p}: {gb.name} exists, kept {len(texts)} records")
        else:
            texts, info = fetch_phylum(p, ncbi_get, a.per_phylum, a.seed, a.pool, a.retmax, acc, gen, a.pause)
            gb.write_text("".join(texts))
            log.append(info)
            print(f"{p}: hits {info['n_hits']}, fetched {info['n_fetched']}, kept {info['n_kept']}", flush=True)
        (a.out / f"refs_{p}.txt").write_text("".join(accession_of(x) + "\n" for x in texts))
    if log:
        f = a.out / "fetch_log.tsv"
        new = not f.exists()
        with open(f, "a", newline="") as h:
            w = csv.DictWriter(h, fieldnames=list(log[0]), delimiter="\t")
            if new:
                w.writeheader()
            w.writerows(log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
