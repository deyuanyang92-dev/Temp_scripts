import csv,sys,collections
from Bio import SeqIO
R="/mnt/n/yzz-分类工作流/projects/mito_benchmark"
rows=[r for r in csv.DictReader(open("/mnt/n/yzz-分类工作流/temp_scripts/Mitoz-annotate/benchmark/sets/mito_bench_2026-10.tsv"),delimiter='\t') if r['split']=='dev']
print('dev records',len(rows))
out=[]
def spans(f):
    return [(int(p.start)+1,int(p.end)) for p in f.location.parts]
for r in rows:
    acc=r['accession']
    try: rec=SeqIO.read(f"{R}/gb/{acc}.gb","genbank")
    except Exception as e: print('skip',acc,e); continue
    L=len(rec.seq)
    feats=[]
    for f in rec.features:
        if f.type in('CDS','tRNA','rRNA'):
            nm=(f.qualifiers.get('gene') or f.qualifiers.get('product') or ['?'])[0]
            sp=spans(f); feats.append((min(s for s,e in sp),max(e for s,e in sp),f.location.strand,f.type,nm,f))
    feats.sort()
    for i,(s,e,st,tp,nm,f) in enumerate(feats):
        if tp!='rRNA': continue
        prod=(f.qualifiers.get('product') or [''])[0].lower(); gn=nm.lower()
        which='rrnS' if ('12s' in prod or 'small' in prod or 'rrns' in gn or '12s' in gn or 'rns' in gn or 's-rrna' in gn) else ('rrnL' if ('16s' in prod or 'large' in prod or 'rrnl' in gn or '16s' in gn or 'rnl' in gn or 'l-rrna' in gn) else '?')
        left=feats[i-1] if i>0 else feats[-1]
        right=feats[i+1] if i+1<len(feats) else feats[0]
        gl=s-left[1]-1 if i>0 else s-1+L-left[1]
        gr=right[0]-e-1 if i+1<len(feats) else L-e+right[0]-1
        out.append((acc,r['phylum'],which,st,s,e,left[4],left[3],gl,right[4],right[3],gr))
import pickle;pickle.dump(out,open('truth_flank.pkl','wb'))
print(len(out))
# gap distributions
for w in ('rrnS','rrnL'):
    for side,idx in (('left',8),('right',11)):
        g=[o[idx] for o in out if o[2]==w]
        c=collections.Counter()
        for x in g:
            c['overlap' if x<0 else ('0' if x==0 else ('1-10' if x<=10 else ('11-50' if x<=50 else ('51-300' if x<=300 else '>300'))))]+=1
        print(w,side,len(g),dict(c))
