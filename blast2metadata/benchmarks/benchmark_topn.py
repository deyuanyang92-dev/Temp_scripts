"""Measure parse + TopN wall time and process peak RSS in isolated subprocesses.

python benchmarks/benchmark_topn.py
python benchmarks/benchmark_topn.py --baseline-module /tmp/old_blasttab.py

The optional baseline is a pre-optimization copy of src/blast2metadata/blasttab.py.
Results must be byte-identical; fixtures intentionally avoid corrected ID/tie bugs.
"""
import argparse
import importlib.util
import json
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-module', type=Path)
    parser.add_argument('--worker', choices=['old', 'new'])
    parser.add_argument('--directory', type=Path)
    args = parser.parse_args()
    if not args.worker:
        with tempfile.TemporaryDirectory(prefix='b2m-topn-bench-') as directory:
            root = Path(directory)
            with (root / 'large.tsv').open('w') as f:
                for i in range(300000):
                    q, s = i % 300, (i // 300) % 100
                    score = (i * 7919) % 1000000 + i / 1000000
                    f.write(f'q{q}\ts{s}\ttitle\t99\t600\t6\t0\t1\t600\t1\t600\t1e-50\t{score}\t100\n')
            for _ in range(3):
                for mode in (['old', 'new'] if args.baseline_module else ['new']):
                    command = [sys.executable, __file__, '--worker', mode, '--directory', directory]
                    if args.baseline_module:
                        command += ['--baseline-module', str(args.baseline_module.resolve())]
                    subprocess.run(command, check=True)
            if args.baseline_module:
                assert (root / 'old.json').read_bytes() == (root / 'new.json').read_bytes()
                print(json.dumps({'comparison': 'byte-identical'}))
        return
    from blast2metadata import blasttab
    module = blasttab
    if args.worker == 'old':
        spec = importlib.util.spec_from_file_location('blast2metadata.old_blasttab', args.baseline_module)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    start = time.perf_counter()
    source = args.directory / 'large.tsv'
    rows = module.parse_blast_results_py(source) if args.worker == 'old' else module.iter_blast_results(source)
    result = module.get_topN_hits_py(rows, 10, 'bitscore')
    elapsed = time.perf_counter() - start
    (args.directory / (args.worker + '.json')).write_text(json.dumps(result, sort_keys=True))
    print(json.dumps(dict(mode=args.worker, seconds=round(elapsed, 3),
                         maxrss_mb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
                         rows=len(result))), flush=True)


if __name__ == '__main__':
    main()
