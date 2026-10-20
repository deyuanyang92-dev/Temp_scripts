import pickle,collections
rows=pickle.load(open('panel_rows.pkl','rb'))
fam={'insect':['Hymenoptera','Lepidoptera','Diptera','Coleoptera','Hemiptera','Orthoptera','Blattodea','Odonata'],
 'vertebrate':['Mammalia','Aves','Actinopterygii','Lepidosauria','Amphibia','Testudines']}
def cls(g):
    for k,v in fam.items():
        if g in v: return k
    return g
def gapbin(x):
    return 'ov' if x<0 else ('0' if x==0 else ('1-10' if x<=10 else ('11-50' if x<=50 else '>50')))
out=collections.defaultdict(list)
for r in rows: out[(cls(r['grp']),r['w'])].append(r)
def cons(seqs,k,end):
    res=[]
    for i in range(k):
        col=[ (s[i] if end==5 else s[-k+i]) for s in seqs if len(s)>=k]
        c=collections.Counter(col); b,n=c.most_common(1)[0]; res.append((b,round(n/len(col),2)))
    return res
for (c,w),rs in sorted(out.items()):
    n=len(rs)
    print(f"== {c} {w} n={n} median len={sorted(r['len'] for r in rs)[n//2]}")
    for end,key in ((5,'n5'),(3,'n3')):
        cnt=collections.Counter((r[key][0],r[key][1].split('_')[0][:8],gapbin(r[key][2])) for r in rs)
        tg=collections.Counter(gapbin(r[key][2]) for r in rs)
        tt=collections.Counter(r[key][0] for r in rs)
        print(f"  {end}' end: neighbour types {dict(tt)} ; gap bins {dict(tg)}")
        print("    top (type,gene,gap):",cnt.most_common(4))
    seqs=[r['seq'] for r in rs]
    print("  5' first10 cons:",''.join(b for b,f in cons(seqs,10,5)),[f for b,f in cons(seqs,10,5)])
    print("  3' last10 cons :",''.join(b for b,f in cons(seqs,10,3)),[f for b,f in cons(seqs,10,3)])
