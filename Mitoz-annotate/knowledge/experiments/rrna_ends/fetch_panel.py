import urllib.request,urllib.parse,json,random,time,sys,csv
random.seed(7)
bench={r['accession'] for r in csv.DictReader(open("/mnt/n/yzz-分类工作流/temp_scripts/Mitoz-annotate/benchmark/sets/mito_bench_2026-10.tsv"),delimiter='\t')}
groups={'Hymenoptera':'Hymenoptera','Lepidoptera':'Lepidoptera','Diptera':'Diptera','Coleoptera':'Coleoptera','Hemiptera':'Hemiptera','Orthoptera':'Orthoptera','Blattodea':'Blattodea','Odonata':'Odonata',
 'Mammalia':'Mammalia','Aves':'Aves','Actinopterygii':'Actinopterygii','Lepidosauria':'Lepidosauria','Amphibia':'Amphibia','Testudines':'Testudines',
 'Chelicerata':'Chelicerata','Crustacea':'Crustacea','Gastropoda':'Gastropoda','Bivalvia':'Bivalvia','Polychaeta':'Polychaeta','Nematoda':'Nematoda','Echinodermata':'Echinodermata','Cnidaria':'Cnidaria','Platyhelminthes':'Platyhelminthes'}
def get(url):
    for i in range(4):
        try: return urllib.request.urlopen(url,timeout=60).read().decode()
        except Exception as e: time.sleep(2)
    raise
for g,org in groups.items():
    term=f'{org}[Organism] AND mitochondrion[Filter] AND refseq[Filter] AND complete genome[Title]'
    u="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&retmax=400&retmode=json&term="+urllib.parse.quote(term)
    ids=json.loads(get(u))['esearchresult']['idlist']; time.sleep(0.4)
    random.shuffle(ids); ids=ids[:16]
    u="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&retmode=text&id="+",".join(ids)
    txt=get(u); open(f'panel/{g}.gb','w').write(txt); print(g,len(ids),len(txt)); time.sleep(0.5)
