import csv,re,sys,collections,statistics
R="/mnt/n/yzz-分类工作流/projects/mito_benchmark/eval/dev/features.tsv"
def parse(s):
    m=re.match(r'(\d+)\.\.(\d+)([+-])$',s)
    if not m: return None
    a,b,st=int(m[1]),int(m[2]),m[3]
    return (a,b,st)
rows=[r for r in csv.DictReader(open(R),delimiter='\t') if r['kind']=='rRNA']
print(collections.Counter((r['tool'],r['class']) for r in rows))
res=collections.defaultdict(list)
for r in rows:
    t=parse(r['truth']); p=parse(r['pred'])
    if not t or not p: 
        res[(r['tool'],r['name'],'unparsed')].append(r); continue
    if t[2]!=p[2]:
        res[(r['tool'],r['name'],'strand')].append(r); continue
    if t[2]=='+':
        d5=t[0]-p[0]; d3=p[1]-t[1]
    else:
        d5=p[1]-t[1]; d3=t[0]-p[0]
    # positive = prediction extends beyond truth (longer)
    res[(r['tool'],r['name'],'ok')].append((r['accession'],d5,d3))
for k,v in sorted(res.items()):
    if k[2]!='ok': print(k,len(v)); continue
    d5=[x[1] for x in v]; d3=[x[2] for x in v]
    ex5=sum(1 for x in d5 if x==0); ex3=sum(1 for x in d3 if x==0); both=sum(1 for x in v if x[1]==0 and x[2]==0)
    w=lambda L,n:sum(1 for x in L if abs(x)<=n)
    print(k,'n',len(v),'5exact',ex5,'3exact',ex3,'both',both,'|5| <=3:',w(d5,3),'3:',w(d3,3),'|<=10:',w(d5,10),w(d3,10),'median d5/d3',statistics.median(d5),statistics.median(d3))
import pickle; pickle.dump(dict(res),open('rr.pkl','wb'))
