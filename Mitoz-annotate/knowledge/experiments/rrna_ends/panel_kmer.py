import pickle,collections
rows=pickle.load(open('panel_rows.pkl','rb'))
def top(seqs,end,k,n=3):
    c=collections.Counter((s[:k] if end==5 else s[-k:]) for s in seqs if len(s)>=k)
    tot=sum(c.values()); return [(m,v,round(v/tot,2)) for m,v in c.most_common(n)], tot
groups=collections.defaultdict(list)
for r in rows: groups[(r['grp'],r['w'])].append(r['seq'])
for (g,w),seqs in sorted(groups.items()):
    if w=='rrnS' or w=='rrnL':
        t5,n=top(seqs,5,5); t3,_=top(seqs,3,5)
        print(g,w,'n',n,"5'5mer",t5,"3'5mer",t3)
