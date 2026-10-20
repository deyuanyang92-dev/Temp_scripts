import json,re,sys,urllib.parse,urllib.request,time
qs=['AFF:"Jiangxi Agricultural University" AND (mitogenome OR "mitochondrial genome" OR "mitochondrial genomes" OR "mitogenomes" OR "mitochondrial DNA")',
    'AFF:"Jiangxi Agricultural University" AND (rRNA OR "ribosomal RNA") AND (mitochondrial)']
seen={}
for q in qs:
  for page in (1,2,3):
    url="https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="+urllib.parse.quote(q)+"&format=json&pageSize=100&page=%d&resultType=lite&sort=P_PDATE_D%%20desc"%page
    d=json.load(urllib.request.urlopen(url,timeout=60))
    for r in d['resultList']['result']:
        seen[r['id']]=r
    if len(d['resultList']['result'])<100: break
print(len(seen))
cand=[r for r in seen.values() if r.get('pmcid')]
print('with pmcid',len(cand))
pat=re.compile(r'[^.]{0,200}(rRNA|rrn[SL]|12S|16S)[^.]{0,200}(boundar|termini|terminus|delimit|secondary structure)[^.]{0,200}\.|[^.]{0,200}(boundar|delimit)[^.]{0,200}(rRNA|rrn[SL]|12S|16S)[^.]{0,200}\.',re.I)
for r in cand:
    try:
        x=urllib.request.urlopen("https://www.ebi.ac.uk/europepmc/webservices/rest/%s/fullTextXML"%r['pmcid'],timeout=60).read().decode()
    except Exception as e: continue
    t=re.sub(r'<[^>]+>',' ',x); t=re.sub(r'\s+',' ',t)
    if not re.search(r'mitochondrial genome|mitogenome',t,re.I): continue
    ms=[m.group(0) for m in pat.finditer(t)]
    if ms:
        print('###',r['pubYear'],r.get('doi'),r['title'][:90],'|',r['authorString'][:40])
        for m in ms[:3]: print('   ',m[:420])
    time.sleep(0.2)
