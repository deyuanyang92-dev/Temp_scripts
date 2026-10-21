"""Fetch per-phylum reference mitogenomes (complete, any database) for rRNA detection experiments.
Up to 150 random records per phylum; benchmark accessions (dev+test, RefSeq and source GenBank) excluded.
Output: refs/<phylum>.gb. Run: python3 -I fetch_refs.py <outdir>"""
import csv, json, random, sys, time, urllib.parse, urllib.request
from pathlib import Path
random.seed(7)
SET = "/mnt/n/yzz-分类工作流/temp_scripts/Mitoz-annotate/benchmark/sets/mito_bench_2026-10.tsv"
rows = list(csv.DictReader(open(SET), delimiter="\t"))
bench = {r["accession"].split(".")[0] for r in rows} | {r["source_accession"].split(".")[0] for r in rows}
phyla = sorted({r["phylum"] for r in rows if r["split"] == "dev"})
out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u):
    for i in range(5):
        try: return urllib.request.urlopen(u, timeout=90).read().decode()
        except Exception: time.sleep(3 * (i + 1))
    raise RuntimeError(u)
for p in phyla:
    f = out / f"{p}.gb"
    if f.exists() and f.stat().st_size: continue
    term = f"{p}[Organism] AND mitochondrion[Filter] AND (complete genome[Title] OR complete sequence[Title]) AND 5000:60000[SLEN]"
    js = json.loads(get(E + "esearch.fcgi?db=nuccore&retmax=5000&retmode=json&term=" + urllib.parse.quote(term)))
    ids = js["esearchresult"]["idlist"]; time.sleep(0.4)
    random.shuffle(ids); ids = ids[:400]
    txt = []
    for i in range(0, len(ids), 100):
        txt.append(get(E + "efetch.fcgi?db=nuccore&rettype=gbwithparts&retmode=text&id=" + ",".join(ids[i:i + 100]))); time.sleep(0.5)
    recs = [r for r in "".join(txt).split("//\n") if r.strip()]
    keep = [r for r in recs if r.split()[1].split(".")[0] not in bench][:150]
    f.write_text("//\n".join(keep) + "//\n")
    print(p, js["esearchresult"]["count"], len(recs), len(keep), flush=True)
