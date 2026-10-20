import pickle,collections
truth=pickle.load(open('truth_flank.pkl','rb'))
def hit(acc,m,mode='g'):
    hs=[]
    for l in open(f'cm/out/{acc}.{m}.{mode}.tbl'):
        if l.startswith('#'): continue
        p=l.split(); hs.append((float(p[15]),min(int(p[7]),int(p[8])),max(int(p[7]),int(p[8])),p[9]))
    return min(hs) if hs else None
def gb(g): return 'ov' if g<0 else ('0' if g==0 else ('1-10' if g<=10 else ('11-50' if g<=50 else '>50')))
for mode in ('g','l'):
    tab=collections.defaultdict(lambda:[0,0,0])
    for (acc,ph,w,st,s,e,ln,lt,gl,rn,rt,gr) in truth:
        if w=='?': continue
        m='mt12' if w=='rrnS' else 'mt16'
        h=hit(acc,m,mode)
        if not h: continue
        _,hs,he,hst=h
        if hst!=('+' if st==1 else '-'): continue
        for side,gap,d in (('L',gl,hs-s),('R',gr,he-e)):
            k=gb(gap); tab[k][0]+=1; tab[k][1]+=(abs(d)<=3); tab[k][2]+=(d==0)
    print('mode',mode,'(g=glocal, l=local) MITOS2 models; truth gap = distance truth end to nearest annotated feature')
    for k in ('ov','0','1-10','11-50','>50'):
        n,w3,ex=tab[k]; print(f"  truth gap {k:6s}: n={n:3d}  CM end within3={w3:3d} ({w3/max(n,1):.2f}) exact={ex}")
