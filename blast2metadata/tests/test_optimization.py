"""Correctness and cache regressions; pure-Python tests also run without BLAST+."""
import argparse
import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from blast2metadata import blasttab, engines, fields, writers


def hit(subject, **kw):
    return dict(qseqid='q', sseqid=subject, bitscore=100., evalue=1e-10,
                pident=99., qcovs=90., **kw)


@pytest.mark.parametrize('threads,jobs,per_job', [(10, 4, None), (7, 4, 5), (3, 2, 12)])
def test_thread_ceiling(threads, jobs, per_job):
    args = argparse.Namespace(threads=threads, parallel_jobs=jobs,
                              local_threads_per_job=per_job, engine='blast')
    pj, tpj = engines.compute_parallel_params(args, 'local', 10)
    assert 1 <= pj * tpj <= threads


def test_distinct_subject_versions_and_local_dotted_ids():
    rows = [hit(s) for s in ['NC_000001.1', 'NC_000001.2', 'contig.partA', 'contig.partB']]
    assert len(blasttab.get_topN_hits_py(rows, 10, 'bitscore')) == 4


def test_ties_independent_of_batch_completion_order():
    rows = [hit('B'), hit('A')]
    assert blasttab.get_topN_hits_py(rows, 1, 'bitscore') == blasttab.get_topN_hits_py(rows[::-1], 1, 'bitscore')


def test_filter_only_requires_enabled_column():
    row = hit('A')
    del row['qcovs']
    assert blasttab.get_topN_hits_py([row], 1, 'bitscore', min_pident=90) == [row]


def test_non_numeric_sort_is_rejected():
    with pytest.raises(ValueError, match='sort_by'):
        blasttab.parse_sort_spec('stitle')


def test_tsv_roundtrip_quoted_titles(tmp_path):
    row = dict(zip(fields.BLAST_HEADERS, ['q', 'A', 'title\twith "quotes"\nand newline', 99., 100, 1, 0, 1, 100, 1, 100, 1e-10, 100., 90.]))
    path = tmp_path / 'hits.tsv'
    writers.rows_to_tsv(path, [row])
    assert blasttab.parse_blast_results_py(path, quoted=True) == [row]


def test_raw_blast_title_with_quotes(tmp_path):
    # BLAST outfmt 6 does not escape quotes; NCBI titles often start with them.
    path = tmp_path / 'raw.tsv'
    path.write_text('q\tA.1\t"Candidatus X" bacterium\t99\t100\t1\t0\t1\t100\t1\t100\t1e-10\t100\t90\n'
                    'q\tB.1\t"unbalanced title\t98\t100\t1\t0\t1\t100\t1\t100\t1e-10\t90\t90\n')
    rows = blasttab.parse_blast_results_py(path)
    assert [r['stitle'] for r in rows] == ['"Candidatus X" bacterium', '"unbalanced title']


def test_bad_table_is_not_silently_empty(tmp_path):
    path = tmp_path / 'broken.tsv'
    path.write_text('q\tA\ttruncated\n')
    with pytest.raises(ValueError, match='line|row|column'):
        blasttab.parse_blast_results_py(path)


def test_standard_12_column_input(tmp_path):
    path = tmp_path / 'std.tsv'
    path.write_text('q\tA\t99\t100\t1\t0\t1\t100\t1\t100\t1e-10\t100\n')
    rows = blasttab.parse_blast_results_py(path)
    assert len(rows) == 1 and rows[0]['qseqid'] == 'q'
    assert blasttab.get_topN_hits_py(rows, 1, 'bitscore')[0]['sseqid'] == 'A'
    with pytest.raises(ValueError, match='qcovs'):
        blasttab.get_topN_hits_py(rows, 1, 'coverage')


def test_streaming_topn_matches_full_sort():
    # Independent full-sort oracle: unsorted input, repeated HSPs, evicted subjects that return.
    rng = random.Random(42)
    rows = [dict(qseqid=f'q{rng.randrange(6)}', sseqid=f's{rng.randrange(80)}',
                 bitscore=rng.random() * 1000, evalue=rng.random(), pident=99., qcovs=90.)
            for _ in range(10000)]
    expected, counts, seen = [], {}, set()
    for r in sorted(rows, key=lambda r: (r['qseqid'], -r['bitscore'], r['evalue'], -r['pident'])):
        pair = r['qseqid'], r['sseqid']
        if pair not in seen and counts.get(pair[0], 0) < 5:
            expected.append(r)
            seen.add(pair)
            counts[pair[0]] = counts.get(pair[0], 0) + 1
    assert blasttab.get_topN_hits_py(iter(rows), 5, 'bitscore') == expected


def test_zero_hit_batch_is_reused(tmp_path, monkeypatch):
    from blast2metadata.cli import parse_arguments
    from blast2metadata.fastautil import BatchInfo
    path = tmp_path / 'q.fa'
    path.write_text('>q\nACGT\n')
    args = parse_arguments(['-i', str(path), '--no_smk', '-bm', 'local', '--resume'])
    calls = []
    monkeypatch.setattr(engines, 'which', lambda _: '/fake/blastn')
    def command(cmd, **kw):
        calls.append(cmd)
        Path(cmd[cmd.index('-out') + 1]).touch()
        return 0, '', ''
    monkeypatch.setattr(engines, 'run_cmd', command)
    batch = BatchInfo('batch_00000', path)
    for _ in range(2):
        assert engines.run_local_blast_one_batch(batch, 'nucleotide', args, tmp_path / 'out', 1, fields.BLAST_HEADERS)[1]
    assert len(calls) == 1


def test_timeout_does_not_publish_partial_output(tmp_path):
    from blast2metadata.util import run_cmd
    output = tmp_path / 'result.part'
    output.write_text('leftover')
    code = 'import pathlib,sys,time; pathlib.Path(sys.argv[1]).write_text("partial"); time.sleep(2)'
    rc, _, err = run_cmd([sys.executable, '-c', code, str(output)], timeout=0.1, output_path=output)
    assert rc == 124 and 'timed out' in err and not output.exists()


def test_cache_rejects_changed_input_parameters_and_modified_output(tmp_path):
    from blast2metadata import cache
    from blast2metadata.cli import parse_arguments
    query = tmp_path / 'q.fa'
    query.write_text('>q\nACTG\n')
    output = tmp_path / 'result.tsv'
    output.touch()  # 零命中是有效的成功输出，但必须带签名。
    args = parse_arguments(['-i', str(query), '--no_smk'])
    sig = cache.search_signature(args, query)
    assert not cache.matches(output, sig)
    cache.commit(output, sig)
    assert cache.matches(output, sig)
    args.evalue = 1e-30
    assert not cache.matches(output, cache.search_signature(args, query))
    query.write_text('>q\nAAAA\n')
    assert not cache.matches(output, cache.search_signature(args, query))
    output.write_text('corrupted')
    assert not cache.matches(output, sig)


def test_translated_blast_selects_subject_database():
    from blast2metadata.cli import parse_arguments
    from blast2metadata.fastautil import get_blast_db
    args = parse_arguments(['-i', 'q.fa', '-bm', 'local', '--blast_exe', '/usr/bin/blastx'])
    assert get_blast_db('nucleotide', args) == 'nr'
    args.blast_exe = '/usr/bin/tblastn'
    assert get_blast_db('protein', args) == 'nt'


def test_metadata_partial_failure_remains_resumable(tmp_path, monkeypatch):
    from blast2metadata import pipeline
    from blast2metadata.cli import parse_arguments
    path = tmp_path / 'ids.txt'
    path.write_text('OK1\nOK2\n')
    args = parse_arguments(['-i', str(path), '-o', str(tmp_path / 'out'), '--mode', 'metadata-only', '--no_smk', '--resume'])
    calls = []
    def fetch(ids, *a, **kw):
        calls.append(ids)
        return {ids[0]: dict(accession=ids[0], organism='test')}
    monkeypatch.setattr(pipeline, 'fetch_metadata_batch', fetch)
    dirs = pipeline.setup_output_directory(Path(args.output))
    assert pipeline.process_one_fasta(str(path), args, dirs) is False
    assert not (Path(args.output) / 'ids_completed.txt').exists()
    assert pipeline.process_one_fasta(str(path), args, dirs) is True
    assert calls == [['OK1', 'OK2'], ['OK2']]


@pytest.mark.parametrize('bad', ['not-a-number', 'nan', 'inf', '1.5'])
def test_corrupt_integer_columns_fail_with_line_number(tmp_path, bad):
    table = tmp_path / 'bad.tsv'
    table.write_text(f'q\tA\t99\t{bad}\t1\t0\t1\t100\t1\t100\t1e-10\t100\n')
    with pytest.raises(ValueError, match='line 1.*length'):
        blasttab.parse_blast_results_py(table)


def test_local_dotted_metadata_id_preserved():
    from blast2metadata.accession import accession_base
    assert accession_base('contig.partA') == 'contig.partA'
    assert accession_base('NC_000001.2') == 'NC_000001'


def test_extra_parameters_cannot_bypass_thread_budget(tmp_path):
    from blast2metadata.cli import parse_arguments
    from blast2metadata.pipeline import preflight_checks
    path = tmp_path / 'q.fa'
    path.write_text('>q\nACTG\n')
    args = parse_arguments(['-i', str(path), '--no_smk', '-bm', 'local', '--extra_para2_blast=-num_threads 100'])
    with pytest.raises(SystemExit, match='Preflight'):
        preflight_checks(args, [str(path)], tmp_path / 'out')
    assert 'managed parameters' in (tmp_path / 'out' / 'preflight_report.json').read_text()


def test_filter_applies_to_best_hsp_not_a_fragment():
    # 主 HSP 85% 未达阈值；同一 subject 的 40 bp/100% 短片段不能让它冒充高一致率命中。
    row = lambda s, bits, pid, length: dict(qseqid='q', sseqid=s, bitscore=bits, evalue=1e-10,
                                            pident=pid, qcovs=95., length=length)
    rows = [row('FALSE.1', 900., 85., 900), row('FALSE.1', 74., 100., 40),
            row('TRUE.1', 500., 98., 300)]
    top = blasttab.get_topN_hits_py(rows, 10, 'bitscore', min_pident=97)
    assert [r['sseqid'] for r in top] == ['TRUE.1']
    # 不设阈值时代表 HSP 仍是 bitscore 最高的那段。
    top = blasttab.get_topN_hits_py(rows, 10, 'bitscore')
    assert [(r['sseqid'], r['length']) for r in top] == [('FALSE.1', 900), ('TRUE.1', 300)]


@pytest.mark.parametrize('extra,expected', [('', 'on'), ('-megablast off', None), ('-megablast on', 'on')])
def test_qblast_blastn_defaults_to_megablast(tmp_path, monkeypatch, extra, expected):
    # NCBI URL API 默认普通 blastn；与 -remote/本地一致应默认 megablast，且可显式关闭。
    from blast2metadata import cli, ratelimit
    from blast2metadata.fastautil import BatchInfo
    q = tmp_path / 'q.fa'
    q.write_text('>q1\nACGTACGTAC\n')
    seen = {}

    def fake_qblast(**kw):
        seen.update(kw)
        raise ConnectionError('stop after capturing parameters')

    monkeypatch.setattr(engines.NCBIWWW, 'qblast', fake_qblast)
    monkeypatch.setattr(ratelimit, 'BLAST_LIMITER', ratelimit._RateLimiter(1000, 1.0))
    args = cli.parse_arguments(['-i', str(q), '--no_smk', '-bm', 'online', '--online_engine', 'biopython',
                                '--retries', '1', '--extra_para2_blast', extra])
    _, ok, _, _ = engines.run_biopython_qblast_one_batch(BatchInfo('b', q), 'nucleotide', args, tmp_path,
                                                         fields.BLAST_HEADERS)
    assert not ok and seen.get('megablast') == expected
