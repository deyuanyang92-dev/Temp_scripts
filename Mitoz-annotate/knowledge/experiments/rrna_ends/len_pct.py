import pickle,collections
rows=pickle.load(open('panel_rows.pkl','rb'))
fam={'insect':['Hymenoptera','Lepidoptera','Diptera','Coleoptera','Hemiptera','Orthoptera','Blattodea','Odonata'],'vertebrate':['Mammalia','Aves','Actinopterygii','Lepidosauria','Amphibia','Testudines']}
def cls(g):
    for k,v in fam.items():
        if g in v: return k
    return g
d=collections.defaultdict(list)
for r in rows: d[(cls(r['grp']),r['w'])].append(r['len'])
import statistics
print("RefSeq panel rRNA lengths (nt) by clade and gene: n, min, 5th pct, median, 95th pct, max")
for k,v in sorted(d.items()):
    v=sorted(v); n=len(v); q=lambda p:v[min(n-1,int(p*n))]
    print(k[0],k[1],n,v[0],q(0.05),v[n//2],q(0.95),v[-1])
