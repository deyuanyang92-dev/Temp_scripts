import sys,collections,statistics
def parse(fn):
    seqs={}
    for l in open(fn):
        if l.startswith('#') or l.startswith('//') or not l.strip(): continue
        n,s=l.split(None,1); seqs[n]=seqs.get(n,'')+s.strip()
    return seqs
def run(fn,label):
    seqs=parse(fn)
    any_=next(iter(seqs.values())); ncol=len(any_)
    midx=[];k=0
    for i in range(ncol):
        c=any_[i]
        if c.isupper() or c=='-': k+=1; midx.append(k)
        else: midx.append(0)
    clen=k
    grp=collections.defaultdict(lambda:[[],[]])
    for n,s in seqs.items():
        cls=n.split('|')[1]
        cols=[midx[i] for i,c in enumerate(s) if c.isupper()]
        if not cols: continue
        grp[cls][0].append(cols[0]); grp[cls][1].append(cols[-1])
    print(f"== {label} CLEN={clen}")
    for cls,(a,b) in sorted(grp.items()):
        n=len(a)
        print(f"  {cls:16s} n={n:3d} start col<=3: {sum(1 for x in a if x<=3):3d} ({sum(1 for x in a if x<=3)/n:.2f}) med {statistics.median(a):6.1f} | end col>=CLEN-3: {sum(1 for x in b if x>=clen-3):3d} ({sum(1 for x in b if x>=clen-3)/n:.2f}) med {statistics.median(b):7.1f}")
for f,l in (('pal/mt_rrnS.sto','MITOS2 rrnS'),('pal/mz_rrnS.sto','MitoZ 12S'),('pal/mt_rrnL.sto','MITOS2 rrnL')):
    run(f,l)
