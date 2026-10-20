import pickle,collections,os
from Bio import SeqIO
R="/mnt/n/yzz-分类工作流/projects/mito_benchmark"
truth=pickle.load(open('truth_flank.pkl','rb'))
ph=pickle.load(open('dev_phylum.pkl','rb'))
T={(a,w):(s,e,st) for (a,p,w,st,s,e,*_) in truth}
def best(f):
    hs=[]
    if not os.path.exists(f): return None
    for l in open(f):
        if l.startswith('#'): continue
        p=l.split(); hs.append((-float(p[14]),min(int(p[7]),int(p[8])),max(int(p[7]),int(p[8])),p[9]))
    return min(hs) if hs else None
cache={}
def nbs(acc):
    if acc not in cache:
        rec=SeqIO.read(f"{R}/gb/{acc}.gb","genbank")
        cache[acc]=[(int(x.location.start)+1,int(x.location.end),x.type) for x in rec.features if x.type in('tRNA','CDS','rRNA')]
    return cache[acc]
def snap(s,e,acc,w,fill):
    ts,te,_=T[(acc,w)]
    nb=[x for x in nbs(acc) if not (x[0]==ts and x[1]==te)]
    ch=True
    while ch:
        ch=False
        for ns,ne,tp in sorted(nb):
            if ne<s or ns>e: continue
            if ns<=s and ne>=e: continue
            if ns<=s: s=ne+1
            elif ne>=e: e=ns-1
            elif ns-s>=e-ne: e=ns-1
            else: s=ne+1
            ch=True
    if fill:
        L=[ne for ns,ne,tp in nb if ne<s]; Rr=[ns for ns,ne,tp in nb if ns>e]
        if L and 0<s-max(L)-1<=fill: s=max(L)+1
        if Rr and 0<min(Rr)-e-1<=fill: e=min(Rr)-1
    return s,e
out=[]
for w,gl,cl in (('rrnS','mt12','S'),('rrnL','mt16','L')):
    for acc,p in ph.items():
        if (acc,w) not in T or not os.path.exists(f'clade/rrn{cl}.{p}.cm'): continue
        ts,te,tst=T[(acc,w)]
        for lab,f in (('universal',f'cm/out/{acc}.{gl}.g.tbl'),('clade',f'clade/{acc}.{cl}.tbl')):
            b=best(f)
            if not b: out.append((w,acc,p,lab,None)); continue
            _,s,e,st=b
            if st!=('+' if tst==1 else '-'): out.append((w,acc,p,lab,'strand')); continue
            rec={'raw':(s,e)}
            for fl in (30,300): rec[f'snap{fl}']=snap(s,e,acc,w,fl)
            out.append((w,acc,p,lab,rec))
pickle.dump(out,open('clade_out.pkl','wb'))
agg=collections.defaultdict(lambda:collections.Counter())
for w,acc,p,lab,rec in out:
    ts,te,_=T[(acc,w)]
    for k in ('raw','snap30','snap300'):
        c=agg[(w,lab,k)]; c['n']+=1
        if not isinstance(rec,dict): continue
        s,e=rec[k]; c['found']+=1; c['L']+=(s==ts); c['R']+=(e==te); c['both']+=(s==ts and e==te)
        c['L3']+=(abs(s-ts)<=3); c['R3']+=(abs(e-te)<=3); c['both3']+=(abs(s-ts)<=3 and abs(e-te)<=3)
for k,c in sorted(agg.items()):
    print(k,dict(c))
# per phylum raw both3
byp=collections.defaultdict(lambda:collections.Counter())
for w,acc,p,lab,rec in out:
    ts,te,_=T[(acc,w)]
    if not isinstance(rec,dict): byp[(w,p,lab)]['miss']+=1; continue
    s,e=rec['raw']; c=byp[(w,p,lab)]; c['n']+=1; c['both3']+=(abs(s-ts)<=3 and abs(e-te)<=3); c['both']+=(s==ts and e==te)
    s,e=rec['snap30']; c['snap30_both']+=(s==ts and e==te)
for k,c in sorted(byp.items()): print(k,dict(c))
