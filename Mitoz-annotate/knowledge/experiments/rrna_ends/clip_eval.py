import pickle,collections,os,statistics,re
exec(open('clade_eval.py').read().split("out=[]")[0])   # reuse T, nbs, snap, ph
def parse_sto(fn):
    seqs={};rf=''
    for l in open(fn):
        if l.startswith('#=GC RF'): rf+=l.split(None,2)[2].strip()
        if l.startswith('#') or l.startswith('//') or not l.strip(): continue
        n,s=l.split(None,1); seqs[n]=seqs.get(n,'')+s.strip()
    return seqs,rf
def panel_cols(fn):
    seqs,_=parse_sto(fn)
    any_=next(iter(seqs.values())); midx=[];k=0
    for c in any_:
        if c.isupper() or c=='-': k+=1; midx.append(k)
        else: midx.append(0)
    d=collections.defaultdict(lambda:[[],[]])
    for n,s in seqs.items():
        g=n.split('|')[2]; cols=[midx[i] for i,c in enumerate(s) if c.isupper()]
        if cols: d[g][0].append(cols[0]); d[g][1].append(cols[-1])
    return d,k
GM={'Arthropoda':['Hymenoptera','Lepidoptera','Diptera','Coleoptera','Hemiptera','Orthoptera','Blattodea','Odonata','Crustacea','Chelicerata'],
 'Chordata':['Mammalia','Aves','Actinopterygii','Lepidosauria','Amphibia','Testudines'],'Mollusca':['Gastropoda','Bivalvia'],'Annelida':['Polychaeta'],
 'Nematoda':['Nematoda'],'Platyhelminthes':['Platyhelminthes'],'Echinodermata':['Echinodermata'],'Cnidaria':['Cnidaria']}
med={}
for w,fn in (('rrnS','pal/mt_rrnS.sto'),('rrnL','pal/mt_rrnL.sto')):
    d,clen=panel_cols(fn)
    for p,gs in GM.items():
        a=[x for g in gs for x in d[g][0]]; b=[x for g in gs for x in d[g][1]]
        med[(w,p)]=(int(statistics.median(a)),int(statistics.median(b)),clen)
print({k:v[:2] for k,v in med.items()})
def hit_map(acc,suffix):
    tbl=[l.split() for l in open(f'aln/{acc}.{suffix}.tbl') if not l.startswith('#')]
    if not tbl: return None
    b=max(tbl,key=lambda p:float(p[14])); s0,e0=int(b[7]),int(b[8])
    seqs,rf=parse_sto(f'aln/{acc}.{suffix}.sto')
    key=f'{acc}/{s0}-{e0}'
    if key not in seqs: return None
    s=seqs[key]; step=1 if e0>=s0 else -1
    pos=s0; col=0; m=[]   # list of (target_pos, modelcol) for residues
    for i,c in enumerate(s):
        ismatch = rf[i]!='.'
        if ismatch: col+=1
        if c in '-.': continue
        m.append((pos,col if ismatch else None)); pos+=step
    return m,('+' if e0>=s0 else '-')
res=collections.defaultdict(collections.Counter); rows=[]
for w,suf in (('rrnS','S'),('rrnL','L')):
    for acc,p in ph.items():
        if (acc,w) not in T or (w,p) not in med: continue
        ts,te,tst=T[(acc,w)]
        hm=hit_map(acc,suf)
        if not hm: res[(w,'miss')]['n']+=1; continue
        m,st=hm
        if st!=('+' if tst==1 else '-'): continue
        c5,c3,clen=med[(w,p)]
        # gene 5' = first residue in alignment order (m is ordered along gene direction since model direction == gene direction)
        keep=[(pos,col) for pos,col in m if col is None or (c5<=col<=c3)]
        # trim leading/trailing inserts
        while keep and keep[0][1] is None: keep.pop(0)
        while keep and keep[-1][1] is None: keep.pop()
        if not keep: continue
        poss=[x[0] for x in keep]; s,e=min(poss),max(poss)
        s0,e0=min(x[0] for x in m),max(x[0] for x in m)
        for lab,(a,b) in (('universal',(s0,e0)),('clip',(s,e)),('clip+snap30',snap(s,e,acc,w,30)),('univ+snap30',snap(s0,e0,acc,w,30)),('clip+snap300',snap(s,e,acc,w,300)),('univ+snap300',snap(s0,e0,acc,w,300))):
            c=res[(w,lab)]; c['n']+=1; c['both']+=(a==ts and b==te); c['both3']+=(abs(a-ts)<=3 and abs(b-te)<=3); c['L']+=(a==ts); c['R']+=(b==te)
for k,c in sorted(res.items()): print(k,dict(c))
