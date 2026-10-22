"""Resume provenance: commit a cache record only after an artifact is complete."""
from __future__ import annotations

import hashlib
import json
import os
import shutil
from pathlib import Path


# 修改 (2026-10-22): 用输入内容 + 语义参数验证缓存，不能仅凭“文件非空”跳过任务。
# 大型数据库只记录文件大小/mtime，不逐字节散列 TB 级库；更换库时须保留正常时间戳。
CACHE_VERSION = 1
SEARCH_OPTIONS = (
    'blast_method', 'engine', 'online_engine', 'blast_exe', 'diamond_exe', 'mmseqs_exe',
    'blast_db', 'blast_db_nt', 'blast_db_nr', 'evalue', 'max_target_seqs',
    'outfmt_blast', 'extra_para2_blast', 'query_gencode', 'batch_size', 'no_blast_header',
)
POST_OPTIONS = (
    'mode', 'topN', 'sort_by', 'min_pident', 'min_qcovs', 'header_style',
    'metadata_type', 'shorter_metadata', 'no_ncbi_metadata', 'local_meta_tsv',
    'local_meta_key', 'local_meta_mode', 'local_meta_prefix', 'pseudo_acc_prefixes',
)


def digest(payload):
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


def file_digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def stamp(path):
    p = Path(path)
    st = p.stat()
    return [str(p.resolve()), st.st_size, st.st_mtime_ns]


def database_stamp(name):
    if not name or name == 'auto':
        return []
    p = Path(name).expanduser()
    roots = [p]
    if not p.is_absolute() and p.parent == Path('.'):
        roots += [Path(d) / p for d in os.environ.get('BLASTDB', '').split(os.pathsep) if d]
    for root in roots:
        # 包括多卷索引/序列文件；无需读取数据库内容。
        files = ([root] if root.is_file() else []) + sorted(root.parent.glob(root.name + '.*'))
        files = [f for f in files if f.is_file()]
        if files:
            return [stamp(f) for f in files]
    return [str(p), 'unresolved']


def search_signature(args, input_path):
    payload = {k: getattr(args, k, None) for k in SEARCH_OPTIONS}
    payload.update(version=CACHE_VERSION, input=file_digest(input_path),
                   input_path=str(Path(input_path).resolve()))
    if getattr(args, 'blast_method', '') == 'local':
        payload['databases'] = {k: database_stamp(getattr(args, k, '')) for k in
                                ('blast_db', 'blast_db_nt', 'blast_db_nr')}
    executables = {}
    for k in ('blast_exe', 'diamond_exe', 'mmseqs_exe'):
        name = getattr(args, k, '')
        for candidate in (['blastn', 'blastp'] if name == 'auto' else [name]):
            exe = shutil.which(candidate) if candidate else None
            if exe:
                executables[candidate] = stamp(exe)
    payload['executables'] = executables
    return digest(payload)


def post_signature(args, search):
    payload = {k: getattr(args, k, None) for k in POST_OPTIONS}
    local = getattr(args, 'local_meta_tsv', '')
    payload.update(search=search, local_meta_content=file_digest(local) if local else None)
    return digest(payload)


def record_path(path):
    return Path(str(path) + '.cache.json')


def matches(path, signature):
    try:
        record = json.loads(record_path(path).read_text())
        return record == dict(signature=signature, artifact=stamp(path))
    except (OSError, ValueError, TypeError):
        return False


def commit(path, signature):
    record = record_path(path)
    part = record.with_name(record.name + '.part')
    part.write_text(json.dumps(dict(signature=signature, artifact=stamp(path)), sort_keys=True))
    os.replace(part, record)


def batch_signature(args, batch):
    # 管道已散列整份输入，批次内复用签名，避免每个 worker 重扫大库和原 FASTA。
    signature = getattr(args, '_search_signature', None)
    return digest([signature or search_signature(args, batch.fasta_path), batch.batch_id])
