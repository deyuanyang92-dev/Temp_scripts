import pickle,csv,collections,sys
from Bio import SeqIO
R="/mnt/n/yzz-分类工作流/projects/mito_benchmark"
truth=pickle.load(open('truth_flank.pkl','rb'))
T={(a,w):(s,e,st) for (a,ph,w,st,s,e,*_) in truth}
def hit(acc,m,mode='g'):
    hs=[]
    for l in open(f'cm/out/{acc}.{m}.{mode}.tbl'):
        if l.startswith('#'): continue
        p=l.split(); hs.append((float(p[15]),min(int(p[7]),int(p[8])),max(int(p[7]),int(p[8])),p[9]))
    return min(hs) if hs else None
def neigh(acc):
    out=[]
    try: rec=SeqIO.read(f"{R}/runs/ours/{acc}/final.gbf","genbank")
    except Exception: return None
    for f in rec.features:
        if f.type in('tRNA','CDS'):
            out.append((int(f.location.start)+1,int(f.location.end),f.type))
    return out
accs=open('cm/dev_acc.txt').read().split()
def sim(model,w,other,fill_max,tonly,mode='g'):
    res=collections.Counter()
    for a in accs:
        if (a,w) not in T: continue
        h=hit(a,model,mode); nb=neigh(a)
        if h is None or nb is None: res['nohit']+=1; continue
        _,s,e,st=h; ts,te,tst=T[(a,w)]
        if st!=('+' if tst==1 else '-'): res['strand']+=1; continue
        oh=hit(a,other,mode); nbs=list(nb)
        if oh: nbs.append((oh[1],oh[2],'rRNA'))
        # trim overlaps
        changed=True
        while changed:
            changed=False
            for ns,ne,tp in sorted(nbs):
                if ne<s or ns>e: continue
                if ns<=s and ne>=e: continue
                if tp=='CDS' and False: continue
                if ns<=s: s=ne+1
                elif ne>=e: e=ns-1
                elif ns-s>=e-ne: e=ns-1
                else: s=ne+1
                changed=True
        if fill_max>0:
            L=[(ne,tp) for ns,ne,tp in nbs if ne<s]; Rr=[(ns,tp) for ns,ne,tp in nbs if ns>e]
            if L:
                ne,tp=max(L)
                if 0<s-ne-1<=fill_max and (not tonly or tp=='tRNA'): s=ne+1
            if Rr:
                ns,tp=min(Rr)
                if 0<ns-e-1<=fill_max and (not tonly or tp=='tRNA'): e=ns-1
        res['n']+=1; res['L']+=(s==ts); res['R']+=(e==te); res['both']+=(s==ts and e==te)
        res['L3']+=(abs(s-ts)<=3); res['R3']+=(abs(e-te)<=3); res['both3']+=(abs(s-ts)<=3 and abs(e-te)<=3)
    return res
for model,w,other in (('mt12','rrnS','mt16'),('mz12','rrnS','mz16'),('mt16','rrnL','mt12'),('mz16','rrnL','mz12')):
    for fm,to in ((0,False),(10,False),(30,False),(50,False),(100,False),(300,False),(300,True),(50,True),(1000,True)):
        r=sim(model,w,other,fm,to)
        print(f"{model} {w} fill<={fm} tRNAonly={to}: n={r['n']} nohit={r['nohit']} L={r['L']} R={r['R']} both={r['both']} | within3 L={r['L3']} R={r['R3']} both={r['both3']}")
