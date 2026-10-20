import pickle,collections
truth=pickle.load(open('truth_flank.pkl','rb'))
def hit(acc,m):
    hs=[]
    for l in open(f'cm/out/{acc}.{m}.g.tbl'):
        if l.startswith('#'): continue
        p=l.split(); hs.append((float(p[15]),min(int(p[7]),int(p[8])),max(int(p[7]),int(p[8])),p[9]))
    return min(hs) if hs else None
tab=collections.defaultdict(lambda:[0,0,0])
for (acc,ph,w,st,s,e,ln,lt,gl,rn,rt,gr) in truth:
    if w=='?': continue
    a=hit(acc,'mt12' if w=='rrnS' else 'mt16'); b=hit(acc,'mz12' if w=='rrnS' else 'mz16')
    if not a or not b: continue
    if a[3]!=b[3] or a[3]!=('+' if st==1 else '-'): continue
    for side,pa,pb,pt in (('L',a[1],b[1],s),('R',a[2],b[2],e)):
        k=(w,'agree<=3' if abs(pa-pb)<=3 else 'disagree')
        tab[k][0]+=1; tab[k][1]+=(abs(pa-pt)<=3); tab[k][2]+=(abs(pb-pt)<=3)
print("MITOS2-CM vs MitoZ-CM raw glocal ends: when the two models' ends agree within 3 nt, how often is MITOS2-CM end within 3 nt of NCBI (col 3) / MitoZ-CM (col 4)")
for k,v in sorted(tab.items()): print(k,'n',v[0],'MITOS2 within3',v[1],f'({v[1]/v[0]:.2f})','MitoZ within3',v[2],f'({v[2]/v[0]:.2f})')
