import sys,subprocess,collections
gene,sto=sys.argv[1],sys.argv[2]   # e.g. rrnS pal/mt_rrnS.sto
C="/root/miniconda3/envs/mitos/bin/"
names=[]
for l in open(sto):
    if l.startswith('#') or l.startswith('//') or not l.strip(): continue
    n=l.split()[0]
    if n not in names: names.append(n)
grp={n:n.split('|')[2] for n in names}
phy={'Nematoda':['Nematoda'],'Platyhelminthes':['Platyhelminthes'],'Echinodermata':['Echinodermata'],'Cnidaria':['Cnidaria'],
 'Annelida':['Polychaeta'],'Mollusca':['Gastropoda','Bivalvia'],
 'Arthropoda':['Hymenoptera','Lepidoptera','Diptera','Coleoptera','Hemiptera','Orthoptera','Blattodea','Odonata','Crustacea','Chelicerata'],
 'Chordata':['Mammalia','Aves','Actinopterygii','Lepidosauria','Amphibia','Testudines']}
for ph,gs in phy.items():
    sel=[n for n in names if grp[n] in gs]
    open(f'clade/{gene}.{ph}.list','w').write('\n'.join(sel)+'\n')
    out=f'clade/{gene}.{ph}.sto'
    with open(out,'w') as fh:
        subprocess.run([C+'esl-alimanip','--seq-k',f'clade/{gene}.{ph}.list','--informat','pfam','--outformat','pfam',sto],stdout=fh,check=True)
    r=subprocess.run([C+'cmbuild','-F',f'clade/{gene}.{ph}.cm',out],capture_output=True,text=True)
    print(ph,len(sel),r.stdout.strip().splitlines()[-3:] if r.returncode==0 else r.stderr[:300])
