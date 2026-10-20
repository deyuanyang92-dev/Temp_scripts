import pickle,collections
from Bio import SeqIO
R="/mnt/n/yzz-分类工作流/projects/mito_benchmark"
truth=pickle.load(open('truth_flank.pkl','rb'))
def hit(acc,m,mode='g'):
    hs=[]
    for l in open(f'cm/out/{acc}.{m}.{mode}.tbl'):
        if l.startswith('#'): continue
        p=l.split(); hs.append((float(p[15]),min(int(p[7]),int(p[8])),max(int(p[7]),int(p[8])),p[9]))
    return min(hs) if hs else None
feats={}
def nbs(acc):
    if acc in feats: return feats[acc]
    rec=SeqIO.read(f"{R}/gb/{acc}.gb","genbank")
    f=[(int(x.location.start)+1,int(x.location.end),x.type) for x in rec.features if x.type in('tRNA','CDS','rRNA')]
    feats[acc]=f; return f
rows=[]
for (acc,ph,w,st,s,e,ln,lt,gl,rn,rt,gr) in truth:
    if w=='?': continue
    m='mt12' if w=='rrnS' else 'mt16'
    h=hit(acc,m)
    if not h: continue
    _,hs,he,hst=h
    if hst!=('+' if st==1 else '-'): continue
    f=[x for x in nbs(acc) if not (x[0]==s and x[1]==e)]
    L=[x for x in f if x[1]<s]; Rr=[x for x in f if x[0]>e]
    if L:
        x=max(L,key=lambda y:y[1]); b=x[1]; rows.append(('L',(b+1)-hs,s-(b+1),hs-s,x[2],w))
    if Rr:
        x=min(Rr,key=lambda y:y[0]); b=x[0]; rows.append(('R',he-(b-1),(b-1)-e,he-e,x[2],w))
def bin_(d):
    for lo,hi,l in ((-10**9,-101,'gap>100'),(-100,-31,'gap31-100'),(-30,-11,'gap11-30'),(-10,-4,'gap4-10'),(-3,-1,'gap1-3'),(0,0,'abut'),(1,3,'ovl1-3'),(4,10,'ovl4-10'),(11,30,'ovl11-30'),(31,10**9,'ovl>30')):
        if lo<=d<=hi: return l
order=['gap>100','gap31-100','gap11-30','gap4-10','gap1-3','abut','ovl1-3','ovl4-10','ovl11-30','ovl>30']
print("MITOS2 glocal CM end vs nearest non-overlapping truth neighbour (bin = distance of CM end to neighbour boundary; ovl = CM end lies inside neighbour)")
print("bin        n   truth_end==neighbour_boundary   truth_end==CM_end   truth_end within3 of CM_end")
tab=collections.defaultdict(lambda:[0,0,0,0])
for side,d,tg,dc,tp,w in rows:
    t=tab[bin_(d)]; t[0]+=1; t[1]+=(tg==0); t[2]+=(dc==0); t[3]+=(abs(dc)<=3)
for k in order:
    t=tab[k]
    if t[0]: print(f"{k:10s} {t[0]:4d}  {t[1]:4d} ({t[1]/t[0]:.2f})   {t[2]:4d} ({t[2]/t[0]:.2f})   {t[3]:4d} ({t[3]/t[0]:.2f})")
