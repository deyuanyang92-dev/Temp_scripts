import pickle,glob,collections,statistics,os
truth=pickle.load(open('truth_flank.pkl','rb'))
T={}
for (acc,ph,w,st,s,e,ln,lt,gl,rn,rt,gr) in truth:
    T[(acc,w)]=(s,e,'+' if st==1 else '-',ph,ln,lt,gl,rn,rt,gr)
def best(acc,m):
    f=f'cm/out/{acc}.{m}.g.tbl'
    hits=[]
    for l in open(f):
        if l.startswith('#'): continue
        p=l.split()
        hits.append((float(p[15]),int(p[5]),int(p[6]),int(p[7]),int(p[8]),p[9],p[16]))
    return min(hits) if hits else None
models={'mz12':'rrnS','mz16':'rrnL','mt12':'rrnS','mt16':'rrnL'}
res={}
for m,w in models.items():
    rows=[]
    for acc in open('cm/dev_acc.txt').read().split():
        if (acc,w) not in T: continue
        b=best(acc,m)
        if b is None: rows.append((acc,None)); continue
        ev,mf,mt,sf,st_,strand,inc=b
        s,e,tst,ph,*_=T[(acc,w)]
        ps,pe=min(sf,st_),max(sf,st_)
        if strand!=tst: rows.append((acc,'strand')); continue
        if tst=='+': d5=s-ps; d3=pe-e
        else: d5=pe-e; d3=s-ps   # careful
        # for minus: 5' of gene is the high coordinate: truth5=e, pred5=pe; extension = pe-e
        if tst=='-': d5=pe-e; d3=s-ps
        rows.append((acc,(d5,d3,mf,mt,ev,inc,ph)))
    res[m]=rows
    ok=[r for r in rows if isinstance(r[1],tuple)]
    d5=[r[1][0] for r in ok]; d3=[r[1][1] for r in ok]
    n=len(rows)
    ex=lambda L:sum(1 for x in L if x==0)
    w3=lambda L:sum(1 for x in L if abs(x)<=3)
    w10=lambda L:sum(1 for x in L if abs(x)<=10)
    both=sum(1 for r in ok if r[1][0]==0 and r[1][1]==0)
    both3=sum(1 for r in ok if abs(r[1][0])<=3 and abs(r[1][1])<=3)
    print(f"{m} n={n} found={len(ok)} 5'exact={ex(d5)} 3'exact={ex(d3)} both={both} both<=3={both3} |5|<=3:{w3(d5)} |3|<=3:{w3(d3)} <=10: {w10(d5)} {w10(d3)} median d5={statistics.median(d5)} d3={statistics.median(d3)}")
    print('   model cols 1st/last:', collections.Counter((r[1][2]==1, r[1][3]) for r in ok).most_common(4))
pickle.dump(res,open('cm_res.pkl','wb'))
