import glob,csv,collections,pickle,re,io
from Bio import SeqIO
bench={r['accession'] for r in csv.DictReader(open("/mnt/n/yzz-分类工作流/temp_scripts/Mitoz-annotate/benchmark/sets/mito_bench_2026-10.tsv"),delimiter='\t')}
rows=[]
for fn in sorted(glob.glob('panel/*.gb')):
    grp=fn.split('/')[-1][:-3]
    for rec in SeqIO.parse(fn,'genbank'):
        acc=rec.id
        if acc in bench: continue
        try: str(rec.seq[:5])
        except Exception: print('noseq',acc); continue
        L=len(rec.seq)
        fs=[]
        for f in rec.features:
            if f.type not in('CDS','tRNA','rRNA','D-loop','misc_feature','rep_origin'): continue
            if len(f.location.parts)>1 and f.type!='CDS': 
                pass
            s=int(f.location.start)+1; e=int(f.location.end)
            nm=(f.qualifiers.get('gene') or f.qualifiers.get('product') or f.qualifiers.get('note') or ['?'])[0]
            fs.append([s,e,f.location.strand,f.type,nm,f])
        fs.sort(key=lambda x:(x[0],x[1]))
        for i,(s,e,st,tp,nm,f) in enumerate(fs):
            if tp!='rRNA': continue
            prod=(f.qualifiers.get('product') or [''])[0].lower(); gn=nm.lower()
            if any(k in prod+gn for k in ('12s','small','rrns','rns','s-rrna','srrna','rrn12')) and '16s' not in prod+gn: w='rrnS'
            elif any(k in prod+gn for k in ('16s','large','rrnl','rnl','l-rrna','lrrna','rrn16')): w='rrnL'
            else: continue
            # neighbours excluding gene-type duplicates
            others=[x for x in fs if x is not fs[i]]
            def nb(side):
                if side=='left':
                    c=[x for x in others if x[1]<s] or others
                    x=max(c,key=lambda y:y[1]); gap=s-x[1]-1
                else:
                    c=[x for x in others if x[0]>e] or others
                    x=min(c,key=lambda y:y[0]); gap=x[0]-e-1
                return x[3],x[4],gap
            nl,nr=nb('left'),nb('right')
            seq=str(rec.seq[s-1:e]).upper()
            if st==-1:
                from Bio.Seq import Seq
                seq=str(Seq(seq).reverse_complement())
                n5,n3=nr,nl   # 5' end of gene is at the right (high coord)
            else: n5,n3=nl,nr
            rows.append(dict(acc=acc,grp=grp,w=w,strand=st,len=e-s+1,seq=seq,n5=n5,n3=n3,org=rec.annotations.get('organism','')))
pickle.dump(rows,open('panel_rows.pkl','wb'))
print(len(rows), collections.Counter((r['grp'],r['w']) for r in rows).most_common(5))
