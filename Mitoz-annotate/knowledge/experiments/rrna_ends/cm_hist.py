import pickle,collections
res=pickle.load(open('cm_res.pkl','rb'))
for m in ('mt12','mz12','mt16','mz16'):
    ok=[r for r in res[m] if isinstance(r[1],tuple)]
    for i,nm in ((0,"5'"),(1,"3'")):
        c=collections.Counter()
        for r in ok:
            d=r[1][i]
            b=('0' if d==0 else ('+1..3' if 0<d<=3 else ('+4..10' if 0<d<=10 else ('+11..50' if 0<d<=50 else ('>+50' if d>50 else ('-1..-3' if -3<=d<0 else ('-4..-10' if -10<=d<0 else ('-11..-50' if -50<=d<0 else '<-50'))))))))
            c[b]+=1
        order=['<-50','-11..-50','-4..-10','-1..-3','0','+1..3','+4..10','+11..50','>+50']
        print(m,nm,len(ok),' '.join(f"{k}:{c[k]}" for k in order))
