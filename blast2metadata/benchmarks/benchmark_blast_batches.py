"""Small local DB benchmark: 1500 queries, <=4 BLAST threads, 3 repetitions.

Run: python benchmarks/benchmark_blast_batches.py
Measures search processes only, excludes DB building/query-file preparation.
Not representative of full nt/nr; compare outputs separately when changing sensitivity.
"""
import concurrent.futures
import json
import subprocess
import time
import tempfile
from pathlib import Path
from Bio import SeqIO
_temp = tempfile.TemporaryDirectory(prefix='b2m-blast-bench-')
root = Path(_temp.name)
root.mkdir(exist_ok=True)
data = Path(__file__).resolve().parents[1] / 'tests' / 'data'
subprocess.run(['makeblastdb','-in',str(data/'mito_db.fasta'),'-dbtype','nucl','-out',str(root/'db')],check=True,capture_output=True)
seqs=list(SeqIO.parse(data/'query.fasta','fasta'))
for repeat in range(3):
    for size,jobs in [(50,4),(500,4),(500,1),(1500,1)]:
        path=root/f'b{size}_{jobs}'
        path.mkdir(exist_ok=True)
        batches=[]
        for start in range(0,1500,size):
            q=path/f'{start:04d}.fa'
            q.write_text(''.join(f'>q{i}\n{seqs[i%len(seqs)].seq}\n' for i in range(start,min(start+size,1500))))
            batches.append(q)
        jobs=min(jobs,len(batches))
        def work(q):
            subprocess.run(['blastn','-query',str(q),'-db',str(root/'db'),'-out',str(q.with_suffix('.tsv')),'-num_threads',str(4//jobs),'-outfmt','6','-max_target_seqs','100'],check=True,capture_output=True)
        t=time.perf_counter()
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
            list(pool.map(work,batches))
        elapsed = time.perf_counter() - t
        # 不同分批方案必须产出相同 HSP 集；顺序不影响正确性。
        rows = sorted(line for q in batches for line in q.with_suffix('.tsv').read_text().splitlines())
        if repeat == 0 and size == 50:
            expected = rows
        assert rows == expected, 'Batch configuration changed BLAST results'
        print(json.dumps(dict(repeat=repeat,batch=size,jobs=jobs,seconds=round(elapsed,3))),flush=True)
