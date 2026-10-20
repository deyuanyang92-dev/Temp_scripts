import sys,collections,re
def parse(fn):
    seqs={}; ss=None
    for l in open(fn):
        if l.startswith('#=GC SS_cons'): ss=ss+l.split()[-1] if ss else l.split()[-1]
        if l.startswith('#') or l.startswith('//') or not l.strip(): continue
        n,s=l.split(None,1); seqs[n]=seqs.get(n,'')+s.strip()
    return seqs,ss
for nm in ('rrnS','rrnL'):
    seqs,ss=parse(f'tal/{nm}.sto')
    any_=next(iter(seqs.values()))
    # match columns: uppercase or '-'
    ncol=len(any_)
    ismatch=[ (c.isupper() or c=='-') for c in any_]
    # verify all seqs agree
    bad=0
    for s in seqs.values():
        for i,c in enumerate(s):
            if c in '.' and ismatch[i]: bad+=1
    midx=[];k=0
    for i in range(ncol):
        if ismatch[i]: k+=1; midx.append(k)
        else: midx.append(0)
    clen=k
    # SS_cons projected on match columns
    ssm=''.join(ss[i] for i in range(ncol) if ismatch[i]) if ss else ''
    paired=[i+1 for i,c in enumerate(ssm) if c in '<>(){}[]']
    print(nm,'CLEN',clen,'bad',bad,'first paired col',paired[0] if paired else None,'last paired col',paired[-1] if paired else None)
    st=[];en=[]
    for name,s in seqs.items():
        cols=[midx[i] for i,c in enumerate(s) if c.isupper()]
        if not cols: continue
        st.append(cols[0]); en.append(cols[-1])
    cs=collections.Counter(st); ce=collections.Counter(en)
    print('  truth start model col: n',len(st),'top',cs.most_common(8))
    print('  truth end model col  : top',ce.most_common(8))
    print('  start col <=3:',sum(1 for x in st if x<=3),' <=10:',sum(1 for x in st if x<=10),' end >=CLEN-3:',sum(1 for x in en if x>=clen-3),' >=CLEN-10:',sum(1 for x in en if x>=clen-10))
    import statistics
    print('  start col median',statistics.median(st),'end col median',statistics.median(en))
    # deleted leading columns by phylum
    byph=collections.defaultdict(lambda:[[],[]])
    for name,s in seqs.items():
        ph=name.split('|')[1]
        cols=[midx[i] for i,c in enumerate(s) if c.isupper()]
        if cols: byph[ph][0].append(cols[0]); byph[ph][1].append(cols[-1])
    for ph,(a,b) in sorted(byph.items()):
        if len(a)>=6: print('   ',ph,len(a),'start cols median',statistics.median(a),'range',min(a),max(a),'| end cols median',statistics.median(b),'range',min(b),max(b))
