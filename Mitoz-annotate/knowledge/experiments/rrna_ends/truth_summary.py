import pickle,collections,csv
truth=pickle.load(open('truth_flank.pkl','rb'))
rows={r['accession']:r for r in csv.DictReader(open("/mnt/n/yzz-分类工作流/temp_scripts/Mitoz-annotate/benchmark/sets/mito_bench_2026-10.tsv"),delimiter='\t')}
def gb(g): return 'ov' if g<0 else ('0' if g==0 else ('1-50' if g<=50 else ('51-300' if g<=300 else '>300')))
t=[x for x in truth if x[2]!='?']
n=len(t)
print('dev rRNAs',n)
print('both ends gap 0:',sum(1 for x in t if x[8]==0 and x[11]==0))
print('both ends |gap|<=10:',sum(1 for x in t if abs(x[8])<=10 and abs(x[11])<=10))
print('both ends gap<=0 (contiguous or overlapping):',sum(1 for x in t if x[8]<=0 and x[11]<=0))
c=collections.Counter()
for x in t:
    c[(x[7],gb(x[8]))]+=1; c[(x[10],gb(x[11]))]+=1
print('per-end neighbour type x gap class (type, class): count')
for ty in ('tRNA','CDS','rRNA'):
    tot=sum(v for (a,b),v in c.items() if a==ty)
    print(' ',ty,'total',tot,{b:c[(ty,b)] for b in ('0','ov','1-50','51-300','>300')})
ph=collections.defaultdict(lambda:[0,0])
for x in t:
    ph[x[1]][0]+=1; ph[x[1]][1]+=(x[8]==0 and x[11]==0)
print('phylum: rRNAs, both ends gap 0')
for k,v in sorted(ph.items()): print(' ',k,v[0],v[1])
tier=collections.defaultdict(lambda:[0,0])
for x in t:
    r=rows[x[0]]; tier[r['evidence_tier']][0]+=1; tier[r['evidence_tier']][1]+=(x[8]==0 and x[11]==0)
print('evidence tier: rRNAs, both ends gap 0')
for k,v in sorted(tier.items()): print(' ',k,v[0],v[1])
