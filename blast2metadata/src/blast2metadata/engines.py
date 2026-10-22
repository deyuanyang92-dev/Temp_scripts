"""Alignment engines: local BLAST+, blastn -remote, Biopython qblast."""
from __future__ import annotations

import argparse
import os
import shutil
import re
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

from Bio import Entrez, SeqIO  # type: ignore
from Bio.Blast import NCBIWWW, NCBIXML  # type: ignore

from . import cache, ratelimit
from .fastautil import BatchInfo, detect_sequence_type, get_blast_db, get_blast_program
from .fields import BLAST_HEADERS
from .fields import build_outfmt_string
from .util import _extra_has_outfmt, _shlex_split, ensure_dir, input_stem, is_file_complete, logger, run_cmd, which

# ─────────────────────────────────────────────
# BLAST runners
# ─────────────────────────────────────────────
def compute_parallel_params(args: argparse.Namespace, mode: str, n_batches: int = 0) -> Tuple[int, int]:
    """(parallel_jobs, threads_per_job). 修改 (2026-10-22): 批次数少于并行数时把线程分给实际批次
    （原来 1 个批次对 nt 只拿到 16/4=4 线程）。"""
    total_threads = max(1, int(args.threads))
    pj = int(args.parallel_jobs) if args.parallel_jobs is not None else 0

    if mode == "online":
        if pj <= 0:
            pj = max(1, int(args.online_max_concurrent))
        pj = max(1, min(pj, int(args.online_max_concurrent)))
        return pj, 1

    if pj <= 0:
        pj = 1 if getattr(args, "engine", "blast") != "blast" else min(4, total_threads)
    pj = max(1, min(pj, total_threads))
    if n_batches > 0:
        pj = min(pj, n_batches)
    if args.local_threads_per_job and int(args.local_threads_per_job) > 0:
        tpj = min(total_threads, int(args.local_threads_per_job))
        pj = min(pj, max(1, total_threads // tpj))
    else:
        # 向下分配保证 pj × tpj <= 总预算；ceil 在 10线程/4任务时会启动12线程。
        tpj = total_threads // pj
    tpj = max(1, tpj)
    return pj, tpj


def run_local_blast_one_batch(
    batch: BatchInfo,
    seq_type: str,
    args: argparse.Namespace,
    blast_dir: Path,
    threads_per_job: int,
    col_names: List[str],
) -> Tuple[str, bool, Path, List[str]]:
    blast_dir = ensure_dir(blast_dir)
    out_file  = blast_dir / f"{batch.batch_id}_blast.tsv"
    # 修改 (2026-10-10): BLAST 先写 .part，成功后改名。理由: is_file_complete 只看大小>0，
    # 中途被杀留下的半截 TSV 会被 --resume 当成已完成。
    part_file = out_file.with_name(out_file.name + ".part")
    signature = cache.batch_signature(args, batch)
    if args.resume and cache.matches(out_file, signature):
        return (batch.batch_id, True, out_file, col_names)

    program = args.blast_exe if args.blast_exe != "auto" else get_blast_program(seq_type)
    if not which(program):
        raise FileNotFoundError(f"BLAST program not found in PATH: {program}")

    db = get_blast_db(seq_type, args)
    # Fix-3: unpack tuple
    outfmt_str, col_names = build_outfmt_string(args.outfmt_blast)

    if _extra_has_outfmt(args.extra_para2_blast):
        raise SystemExit(
            "ERROR: --extra_para2_blast must NOT include -outfmt. "
            "Remove it and use --outfmt_blast instead."
        )

    cmd = [
        program,
        "-query", str(batch.fasta_path),
        "-db",    str(db),
        "-outfmt", outfmt_str,
        "-out",   str(part_file),
        "-evalue", str(args.evalue),
        "-max_target_seqs", str(args.max_target_seqs),
        "-num_threads", str(max(1, int(threads_per_job))),
    ]
    if Path(program).name in {"blastx", "tblastx"} and getattr(args, "query_gencode", None):
        cmd += ["-query_gencode", str(args.query_gencode)]
    if args.extra_para2_blast:
        cmd += _shlex_split(args.extra_para2_blast)

    rc, so, se = run_cmd(cmd, retries=args.retries, retry_sleep=args.retry_sleep,
                         timeout=getattr(args, "blast_timeout", 0), output_path=part_file)
    if rc != 0:
        logger.error(f"[LOCAL_BLAST_FAIL] {batch.batch_id}\n{se[-2000:]}")
        return (batch.batch_id, False, out_file, col_names)
    part_file.touch()  # BLAST+ may skip -out when there are no hits
    os.replace(part_file, out_file)
    cache.commit(out_file, signature)
    return (batch.batch_id, True, out_file, col_names)


def run_remote_blast_one_batch(
    batch: BatchInfo,
    seq_type: str,
    args: argparse.Namespace,
    blast_dir: Path,
    col_names: List[str],
) -> Tuple[str, bool, Path, List[str]]:
    blast_dir = ensure_dir(blast_dir)
    out_file  = blast_dir / f"{batch.batch_id}_blast.tsv"
    # 修改 (2026-10-10): BLAST 先写 .part，成功后改名。理由: is_file_complete 只看大小>0，
    # 中途被杀留下的半截 TSV 会被 --resume 当成已完成。
    part_file = out_file.with_name(out_file.name + ".part")
    signature = cache.batch_signature(args, batch)
    if args.resume and cache.matches(out_file, signature):
        return (batch.batch_id, True, out_file, col_names)

    program = args.blast_exe if args.blast_exe != "auto" else get_blast_program(seq_type)
    if not which(program):
        raise FileNotFoundError(
            f"BLAST+ program not found in PATH: {program}\n"
            "Install BLAST+ (blastn/blastp) or use --blast_method local."
        )

    db = get_blast_db(seq_type, args)
    outfmt_str, col_names = build_outfmt_string(args.outfmt_blast)

    if _extra_has_outfmt(args.extra_para2_blast):
        raise SystemExit(
            "ERROR: --extra_para2_blast must NOT include -outfmt. "
            "Remove it and use --outfmt_blast instead."
        )

    cmd = [
        program, "-remote",
        "-query", str(batch.fasta_path),
        "-db",    str(db),
        "-outfmt", outfmt_str,
        "-out",   str(part_file),
        "-evalue", str(args.evalue),
        "-max_target_seqs", str(args.max_target_seqs),
    ]
    if args.extra_para2_blast:
        cmd += _shlex_split(args.extra_para2_blast)

    # Fix-G: use BLAST_LIMITER for remote BLAST calls
    ratelimit.BLAST_LIMITER.acquire()
    rc, so, se = run_cmd(cmd, retries=args.retries, retry_sleep=args.retry_sleep,
                         timeout=getattr(args, "blast_timeout", 0), output_path=part_file)
    if rc != 0:
        logger.error(f"[REMOTE_BLAST_FAIL] {batch.batch_id}\n{se[-2000:]}")
        return (batch.batch_id, False, out_file, col_names)
    part_file.touch()
    os.replace(part_file, out_file)
    cache.commit(out_file, signature)

    if args.online_delay_sec > 0:
        time.sleep(float(args.online_delay_sec))
    return (batch.batch_id, True, out_file, col_names)


def _parse_extra_params_for_qblast(extra: str) -> Dict[str, Union[str, bool]]:
    params: Dict[str, Union[str, bool]] = {}
    if not extra:
        return params
    xs = _shlex_split(extra)
    i = 0
    while i < len(xs):
        if xs[i].startswith("-"):
            k = xs[i].lstrip("-")
            if i + 1 < len(xs) and not xs[i + 1].startswith("-"):
                params[k] = xs[i + 1]
                i += 2
            else:
                params[k] = True
                i += 1
        else:
            i += 1
    return params


def entrez_setup(args: argparse.Namespace) -> None:
    # 修改 (2026-10-22): qblast 走 NCBIWWW.email/tool，原代码只设 Entrez.*，qblast 请求不带 email。
    if args.email:
        Entrez.email = args.email
        NCBIWWW.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
    if args.entrez_tool:
        Entrez.tool = args.entrez_tool
        NCBIWWW.tool = args.entrez_tool


def _count_gap_opens(aligned: str) -> int:
    return len(re.findall(r"-+", aligned or ""))


def _qcovs_from_hsps(hsps: Sequence[object], qlen: int) -> int:
    """BLAST+ qcovs: % of query covered by the union of all HSPs to one subject."""
    spans = [
        (min(h.query_start, h.query_end), max(h.query_start, h.query_end))  # type: ignore[attr-defined]
        for h in hsps
    ]
    return _qcovs_from_spans(spans, qlen)


def _qblast_rows(br: object, qid: str, qlen: int) -> List[str]:
    """One BLAST XML record -> outfmt-6-style lines in BLAST_HEADERS order."""
    # 修改 (2026-10-10): 列含义对齐 BLAST+ outfmt 6。
    # 理由: 原代码 bitscore 写 hsp.score（原始分），gapopen 写总 gap 数，
    #       qcovs 写单个 HSP 覆盖度；与本地 BLAST 结果不可比。
    lines = []
    for aln in br.alignments:  # type: ignore[attr-defined]
        subj_id = getattr(aln, "accession", "") or getattr(aln, "hit_id", "")
        title   = getattr(aln, "title", "") or ""
        qcovs   = _qcovs_from_hsps(aln.hsps, qlen)
        for hsp in aln.hsps:
            pid     = (hsp.identities / max(1, hsp.align_length)) * 100.0
            mism    = hsp.align_length - hsp.identities - (getattr(hsp, "gaps", 0) or 0)
            gapopen = _count_gap_opens(hsp.query) + _count_gap_opens(hsp.sbjct)
            lines.append(
                f"{qid}\t{subj_id}\t{title}\t{pid:.3f}\t{hsp.align_length}\t"
                f"{mism}\t{gapopen}\t{hsp.query_start}\t{hsp.query_end}\t"
                f"{hsp.sbjct_start}\t{hsp.sbjct_end}\t"
                f"{hsp.expect}\t{hsp.bits}\t{qcovs}\n"
            )
    return lines


def run_biopython_qblast_one_batch(
    batch: BatchInfo,
    seq_type: str,
    args: argparse.Namespace,
    blast_dir: Path,
    col_names: List[str],
) -> Tuple[str, bool, Path, List[str]]:
    blast_dir = ensure_dir(blast_dir)
    out_file  = blast_dir / f"{batch.batch_id}_blast.tsv"
    signature = cache.batch_signature(args, batch)
    if args.resume and cache.matches(out_file, signature):
        return (batch.batch_id, True, out_file, col_names)

    if _extra_has_outfmt(args.extra_para2_blast):
        raise SystemExit(
            "ERROR: --extra_para2_blast must NOT include -outfmt. "
            "Remove it and use --outfmt_blast instead."
        )

    part_file = out_file.with_name(out_file.name + ".part")
    entrez_setup(args)
    qblast_program = get_blast_program(seq_type) if args.blast_exe == "auto" else Path(args.blast_exe).name
    db    = get_blast_db(seq_type, args)
    extra = _parse_extra_params_for_qblast(args.extra_para2_blast)
    # 修改 (2026-10-23): blastn 默认 MEGABLAST=on。
    # 理由: NCBI URL API 不带 MEGABLAST=on 时跑普通 blastn (word 11, 2/-3 计分)，而 blastn -remote、
    #       本地 BLAST+ 和 NCBI 网页默认都是 megablast，三条路径的 bitscore 不可比。实测 10 条 COI
    #       对 nt: 普通 blastn 125.2 s，megablast 27.7 s；最佳命中 pident 全部相同，但计分不同使
    #       多个同分 subject 的 TopN 成员变化（q5 前 10 只重合 2 个）。需要更敏感的普通 blastn 时
    #       用 --extra_para2_blast "-megablast off"（此处删掉该参数，交给 URL API 默认）。
    if qblast_program == "blastn":
        if str(extra.get("megablast", "on")).lower() in {"off", "false", "no", "0"}:
            extra.pop("megablast", None)
        else:
            extra["megablast"] = "on"

    seqs = list(SeqIO.parse(str(batch.fasta_path), "fasta"))

    # 修改 (2026-10-23): 一个批次只提交一次 qblast（多条 FASTA），不再逐条提交。
    # 理由: Biopython 1.85 的 qblast 提交后先等 20 s 才第一次取结果，之后每 60 s 一次
    #       （NCBI 规定同一 RID 每分钟最多查询一次），再加原来每条之后固定 sleep 12 s，
    #       单条序列至少 30–90 s，1000 条要 9–25 小时。NCBI 的 URL API 支持一次提交多条
    #       query，XML 里每条 query 一个 Iteration，所以等待次数按批而不是按条计算。
    # 原代码: for rec in seqs: NCBIWWW.qblast(sequence=str(rec.seq), ...); time.sleep(online_delay_sec)
    def search(recs: List[object], label: str) -> List[str]:
        fasta = "".join(f">{r.id}\n{r.seq}\n" for r in recs)  # type: ignore[attr-defined]
        want = [r.id for r in recs]  # type: ignore[attr-defined]
        delay = 5
        last: Optional[Exception] = None
        for attempt in range(1, args.retries + 1):
            try:
                logger.info(f"[QBLAST] {batch.batch_id} {label} n={len(recs)} attempt {attempt}")
                # 提交前限速；之后的轮询由 Biopython 自己按 NCBI 规定节流。
                ratelimit.BLAST_LIMITER.acquire()
                handle = NCBIWWW.qblast(
                    program=qblast_program,
                    database=db,
                    sequence=fasta,
                    expect=args.evalue,
                    hitlist_size=args.max_target_seqs,
                    **extra,
                )
                try:
                    records = list(NCBIXML.parse(handle))
                finally:
                    handle.close()
                # 结果按提交顺序返回，每条 query 一个 record（无命中也有）；按 ID 核对，
                # 防止 NCBI 漏掉或重排某条 query 时把命中错配到别的序列上。
                got = [str(br.query).split()[0] if br.query else "" for br in records]
                if got != want:
                    raise ValueError(f"qblast returned queries {got[:3]}... for {want[:3]}... ({len(got)}/{len(want)})")
                return [line for br, r in zip(records, recs) for line in _qblast_rows(br, r.id, len(r.seq))]  # type: ignore[attr-defined]
            except Exception as e:
                last = e
                logger.warning(f"[QBLAST_FAIL] {batch.batch_id} {label} n={len(recs)} attempt {attempt}: {e}")
                if attempt < args.retries:
                    time.sleep(delay)
                    delay = int(delay * 1.8)
        # 整批被 NCBI 拒绝（Biopython 对 NCBI 错误页抛 ValueError，如单次请求 CPU 超限、
        # 某条序列不合法；或上面的 ID 核对失败）时对半拆分重试，只有单条也失败才算失败。
        # 网络错误不拆分：断网时拆分只会把请求数放大到约 2n 倍，留给 --resume 重跑即可。
        if len(recs) > 1 and isinstance(last, ValueError):
            half = len(recs) // 2
            logger.info(f"[QBLAST] {batch.batch_id} {label}: splitting {len(recs)} queries and retrying")
            return search(recs[:half], label + "a") + search(recs[half:], label + "b")
        raise RuntimeError(f"qblast failed for {want[0]}: {last}")

    try:
        lines = search(seqs, "all") if seqs else []
    except RuntimeError as e:
        logger.error(f"[QBLAST_FAIL] {batch.batch_id}: {e}")
        return (batch.batch_id, False, out_file, col_names)
    with open(part_file, "w", encoding="utf-8") as out:
        out.writelines(lines)
    if args.online_delay_sec > 0:
        time.sleep(float(args.online_delay_sec))

    os.replace(part_file, out_file)
    cache.commit(out_file, signature)
    return (batch.batch_id, True, out_file, col_names)


# ─────────────────────────────────────────────
# DIAMOND / MMseqs2 local engines (2026-10-22)
# ─────────────────────────────────────────────
# Both tools write BLAST-tabular output, so the rest of the pipeline (TopN, metadata, headers)
# is unchanged. Differences handled here:
#   • neither tool has BLAST's per-subject `qcovs`; we request `qlen` and compute qcovs from the
#     union of all HSPs per query-subject pair, exactly like BLAST+ (see _qcovs_from_hsps).
#   • DIAMOND only searches protein databases: protein query -> blastp, nucleotide query -> blastx.
#   • MMseqs2 needs --search-type (1 aa, 2 translated, 3 nucleotide) and --strand 2 for DNA/DNA;
#     the target DB type is read from <db>.dbtype.
ENGINE_CHOICES = ["blast", "diamond", "mmseqs"]
_ENGINE_FASTA_EXTS = (".fa", ".fasta", ".fas", ".faa", ".fna", ".pep", ".aa")
_DIAMOND_COLS = BLAST_HEADERS[:-1] + ["qlen"]  # qcovs replaced by qlen, computed afterwards
_MMSEQS_COLS = ("query,target,theader,pident,alnlen,mismatch,gapopen,"
                "qstart,qend,tstart,tend,evalue,bits,qlen")


def _is_fasta_file(path: str) -> bool:
    p = path.lower()
    if p.endswith(".gz"):
        p = p[:-3]
    return os.path.isfile(path) and p.endswith(_ENGINE_FASTA_EXTS)


def _fasta_is_protein(path: str) -> bool:
    return detect_sequence_type(path) == "protein"


def mmseqs_db_type(db_prefix: str) -> Optional[str]:
    """'nucleotide' | 'protein' from <db>.dbtype (lower 16 bits: 0 = amino acid, 1 = nucleotide)."""
    try:
        with open(db_prefix + ".dbtype", "rb") as fh:
            v = int.from_bytes(fh.read(4), "little") & 0xFFFF
    except OSError:
        return None
    return {0: "protein", 1: "nucleotide"}.get(v)


def _db_lock(path: Path):
    """Advisory lock so parallel Snakemake children never build the same DB twice."""
    import fcntl

    class _Lock:
        def __enter__(self):
            ensure_dir(path.parent)
            self.fh = open(path, "w")
            fcntl.flock(self.fh, fcntl.LOCK_EX)
            return self

        def __exit__(self, *a):
            fcntl.flock(self.fh, fcntl.LOCK_UN)
            self.fh.close()

    return _Lock()


_BLASTDB_SUFFIXES = {"protein": (".pal", ".pin", ".00.pin"), "nucleotide": (".nal", ".nin", ".00.nin")}
DIAMOND_MIN_BLASTDB_VERSION = (2, 1, 14)  # wiki: "downloading and using a BLAST database (use DIAMOND >= v2.1.14)"


def find_blast_db(name: str, db_type: str) -> Optional[str]:
    """Absolute base path of a BLAST database (`nr`, `/db/nr`, `nr.pal`), searching $BLASTDB for bare
    names like blastn does; None if not found. DIAMOND itself does not consult $BLASTDB."""
    name = os.path.expanduser(name)
    for suf in (".pal", ".nal", ".pin", ".nin"):
        if name.endswith(suf):
            name = name[: -len(suf)]
    dirs = [""] if (os.sep in name or os.path.isabs(name)) else [""] + [
        d for d in os.environ.get("BLASTDB", "").split(os.pathsep) if d]
    for d in dirs:
        base = os.path.join(d, name) if d else name
        if any(os.path.isfile(base + suf) for suf in _BLASTDB_SUFFIXES[db_type]):
            return os.path.abspath(base)
    return None


def diamond_version(exe: str) -> Tuple[int, ...]:
    rc, so, se = run_cmd([exe, "version"], retries=1)
    m = re.search(r"version\s+(\d+)\.(\d+)\.(\d+)", so + se)
    return tuple(int(x) for x in m.groups()) if m else (0,)


def resolve_engine_db(engine: str, requested: str, seq_type: str, out_dir: Path,
                      exe: str, threads: int = 1) -> Tuple[str, str]:
    """Return (db_path, db_seq_type) for DIAMOND/MMseqs2, building the DB from a FASTA if needed.

    `requested` may be a ready DB (x.dmnd; a BLAST protein DB such as nr for DIAMOND >= 2.1.14;
    an MMseqs2 DB prefix with x.dbtype, e.g. from `mmseqs databases NR`) or a FASTA file; a FASTA is
    converted once into OUT/_engine_db/ (guarded by a lock, so --resume and parallel jobs reuse it).
    """
    req = os.path.expanduser(requested)
    if engine == "diamond":
        for cand in (req, req + ".dmnd"):
            if cand.endswith(".dmnd") and os.path.isfile(cand):
                return os.path.abspath(cand), "protein"
        # 修改 (2026-10-22): DIAMOND >= 2.1.14 直接搜 BLAST (v5) 蛋白库，nr 用
        # `update_blastdb.pl --decompress --blastdb_version 5 nr` 下载后无需 makedb。
        blast_prot = find_blast_db(req, "protein")
        if blast_prot:
            ver = diamond_version(exe)
            if ver < DIAMOND_MIN_BLASTDB_VERSION:
                raise SystemExit(
                    f"ERROR: {exe} is DIAMOND {'.'.join(map(str, ver))}; searching a BLAST database directly "
                    f"needs DIAMOND >= 2.1.14. Upgrade, or build a .dmnd with `diamond makedb`."
                )
            if not any(os.path.isfile(blast_prot + s) for s in (".pdb", ".00.pdb", ".pos")):
                logger.warning(f"[DIAMOND] {blast_prot} looks like a BLAST v4 database; NCBI nr downloads "
                               "(`update_blastdb.pl --blastdb_version 5 nr`) are v5, which DIAMOND documents.")
            return blast_prot, "protein"
        if find_blast_db(req, "nucleotide"):
            raise SystemExit(
                f"ERROR: --engine diamond: {requested} is a nucleotide BLAST database. DIAMOND searches "
                "protein databases only (blastp/blastx); use --engine blast or --engine mmseqs for nt."
            )
        if not _is_fasta_file(req):
            raise SystemExit(
                f"ERROR: --engine diamond needs a protein database: a BLAST protein DB (e.g. nr, DIAMOND >= 2.1.14), "
                f"a .dmnd file, or a protein FASTA; got: {requested}\n"
                "Pass it with --blast_db. DIAMOND searches protein databases only (blastp / blastx)."
            )
        if not _fasta_is_protein(req):
            raise SystemExit(f"ERROR: --engine diamond: {requested} looks like a nucleotide FASTA; "
                             "DIAMOND databases must be protein.")
        db_dir = ensure_dir(out_dir.resolve() / "_engine_db")
        target = db_dir / (input_stem(req) + ".dmnd")
        with _db_lock(db_dir / ".lock"):
            if not is_file_complete(target):
                logger.info(f"[DIAMOND] building database {target} from {req}")
                rc, so, se = run_cmd([exe, "makedb", "--in", req, "-d", str(target.with_suffix("")),
                                      "-p", str(max(1, threads)), "--quiet"], retries=1)
                if rc != 0 or not is_file_complete(target):
                    raise SystemExit(f"ERROR: diamond makedb failed (rc={rc}):\n{se[-1500:]}")
        return str(target), "protein"

    if engine == "mmseqs":
        for cand in (req, req[:-7] if req.endswith(".dbtype") else req):
            t = mmseqs_db_type(cand)
            if t:
                return os.path.abspath(cand), t
        if not _is_fasta_file(req):
            hint = ""
            if find_blast_db(req, "protein") or find_blast_db(req, "nucleotide"):
                hint = (" MMseqs2 cannot read BLAST database files; download its own copy with "
                        "`mmseqs databases NR|NT|core_nt <outDB> <tmp>` or point --blast_db at a FASTA.")
            raise SystemExit(
                f"ERROR: --engine mmseqs needs an MMseqs2 database (`mmseqs createdb` / `mmseqs databases`) "
                f"or a FASTA; got: {requested}. Pass it with --blast_db." + hint
            )
        db_dir = ensure_dir(out_dir.resolve() / "_engine_db")
        target = db_dir / (input_stem(req) + "_mmseqsdb")
        with _db_lock(db_dir / ".lock"):
            if mmseqs_db_type(str(target)) is None:
                logger.info(f"[MMSEQS] building database {target} from {req}")
                rc, so, se = run_cmd([exe, "createdb", req, str(target), "-v", "1"], retries=1)
                t = mmseqs_db_type(str(target))
                if rc != 0 or t is None:
                    raise SystemExit(f"ERROR: mmseqs createdb failed (rc={rc}):\n{se[-1500:]}")
        return str(target), mmseqs_db_type(str(target)) or "nucleotide"

    raise ValueError(engine)


def prepare_engine_db(args: argparse.Namespace, seq_type: str, out_dir: Path) -> None:
    """Resolve/build the DIAMOND or MMseqs2 database once and store it on args (serialised to children)."""
    if args.blast_method != "local" or args.engine == "blast" or getattr(args, "_engine_db", ""):
        return
    if args.engine == "diamond":
        requested = args.blast_db if args.blast_db != "auto" else args.blast_db_nr
        exe = args.diamond_exe
    else:
        requested = get_blast_db(seq_type, args)
        exe = args.mmseqs_exe
    if not which(exe):
        raise SystemExit(f"ERROR: --engine {args.engine}: executable not found: {exe}")
    db, db_type = resolve_engine_db(args.engine, requested, seq_type, out_dir, exe, threads=int(args.threads))
    args._engine_db = db
    args._engine_db_type = db_type
    logger.info(f"[ENGINE] {args.engine} db={db} db_type={db_type} (requested: {requested})")


def mmseqs_search_type(query_type: str, db_type: str) -> int:
    if query_type == "nucleotide" and db_type == "nucleotide":
        return 3
    if query_type == "nucleotide" and db_type == "protein":
        return 2
    if query_type == "protein" and db_type == "protein":
        return 1
    raise SystemExit("ERROR: --engine mmseqs: protein queries against a nucleotide database are not "
                     "supported (tblastn-like search); use BLAST+ tblastn via --engine blast --blast_exe tblastn.")


def _rewrite_with_qcovs(part_file: Path, has_header: bool = False) -> None:
    """Replace the trailing `qlen` column by BLAST-style `qcovs` (union of HSPs per query-subject)."""
    rows: List[List[str]] = []
    with open(part_file, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(_DIAMOND_COLS):
                continue
            # stitle may contain tabs: keep the first 2 and the last 11 columns, join the middle
            if len(parts) > len(_DIAMOND_COLS):
                tail = len(_DIAMOND_COLS) - 3
                parts = parts[:2] + ["\t".join(parts[2:len(parts) - tail])] + parts[len(parts) - tail:]
            rows.append(parts)
    spans: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    qlens: Dict[Tuple[str, str], int] = {}
    for r in rows:
        key = (r[0], r[1])
        a, b = int(float(r[7])), int(float(r[8]))
        spans.setdefault(key, []).append((min(a, b), max(a, b)))
        qlens[key] = int(float(r[13]))
    qcovs = {k: _qcovs_from_spans(v, qlens[k]) for k, v in spans.items()}
    with open(part_file, "w", encoding="utf-8") as out:
        for r in rows:
            out.write("\t".join(r[:13] + [str(qcovs[(r[0], r[1])])]) + "\n")


def _qcovs_from_spans(spans: List[Tuple[int, int]], qlen: int) -> int:
    """BLAST+ qcovs: % of query covered by the union of HSP spans. Rounds half up but reports 100
    only for full coverage (checked against blastn 2.17)."""
    covered, cur_s, cur_e = 0, None, None
    for a, b in sorted(spans):
        if cur_e is None or a > cur_e + 1:
            if cur_e is not None:
                covered += cur_e - cur_s + 1
            cur_s, cur_e = a, b
        else:
            cur_e = max(cur_e, b)
    if cur_e is not None:
        covered += cur_e - cur_s + 1
    pct = int(covered / max(1, qlen) * 100 + 0.5)
    return 99 if pct == 100 and covered < qlen else pct


def run_diamond_one_batch(
    batch: BatchInfo, seq_type: str, args: argparse.Namespace, blast_dir: Path,
    threads_per_job: int, col_names: List[str],
) -> Tuple[str, bool, Path, List[str]]:
    blast_dir = ensure_dir(blast_dir)
    out_file = blast_dir / f"{batch.batch_id}_blast.tsv"
    part_file = out_file.with_name(out_file.name + ".part")
    signature = cache.batch_signature(args, batch)
    if args.resume and cache.matches(out_file, signature):
        return (batch.batch_id, True, out_file, col_names)
    exe = args.diamond_exe
    if not which(exe):
        raise FileNotFoundError(f"DIAMOND executable not found: {exe}")
    program = "blastp" if seq_type == "protein" else "blastx"
    cmd = [
        exe, program,
        "-q", str(batch.fasta_path),
        "-d", str(args._engine_db),
        "-o", str(part_file),
        "-f", "6", *_DIAMOND_COLS,
        "-e", str(args.evalue),
        "-k", str(args.max_target_seqs),
        "-p", str(max(1, int(threads_per_job))),
        "--quiet",
    ]
    if program == "blastx" and getattr(args, "query_gencode", None):
        cmd += ["--query-gencode", str(args.query_gencode)]
    if args.extra_para2_blast:
        cmd += _shlex_split(args.extra_para2_blast)
    rc, so, se = run_cmd(cmd, retries=args.retries, retry_sleep=args.retry_sleep,
                         timeout=getattr(args, "blast_timeout", 0), output_path=part_file)
    if rc != 0:
        logger.error(f"[DIAMOND_FAIL] {batch.batch_id}\n{se[-2000:]}")
        return (batch.batch_id, False, out_file, col_names)
    part_file.touch()
    _rewrite_with_qcovs(part_file)
    os.replace(part_file, out_file)
    cache.commit(out_file, signature)
    return (batch.batch_id, True, out_file, col_names)


def run_mmseqs_one_batch(
    batch: BatchInfo, seq_type: str, args: argparse.Namespace, blast_dir: Path,
    threads_per_job: int, col_names: List[str],
) -> Tuple[str, bool, Path, List[str]]:
    blast_dir = ensure_dir(blast_dir)
    out_file = blast_dir / f"{batch.batch_id}_blast.tsv"
    part_file = out_file.with_name(out_file.name + ".part")
    signature = cache.batch_signature(args, batch)
    if args.resume and cache.matches(out_file, signature):
        return (batch.batch_id, True, out_file, col_names)
    exe = args.mmseqs_exe
    if not which(exe):
        raise FileNotFoundError(f"MMseqs2 executable not found: {exe}")
    search_type = mmseqs_search_type(seq_type, str(args._engine_db_type))
    tmp_dir = batch.fasta_path.parent / f"{batch.batch_id}_mmseqs_tmp"
    k = int(args.max_target_seqs)
    cmd = [
        exe, "easy-search", str(batch.fasta_path), str(args._engine_db), str(part_file), str(tmp_dir),
        "--search-type", str(search_type),
        "--format-output", _MMSEQS_COLS,
        "-e", str(args.evalue),
        "--max-seqs", str(max(300, k)),
        "--max-accept", str(k),
        "--threads", str(max(1, int(threads_per_job))),
        "-v", "1",
    ]
    if search_type == 3:
        cmd += ["--strand", "2"]  # BLAST searches both strands; MMseqs2 defaults to forward only
    if search_type == 2 and getattr(args, "query_gencode", None):
        cmd += ["--translation-table", str(args.query_gencode)]
    if args.extra_para2_blast:
        cmd += _shlex_split(args.extra_para2_blast)
    rc, so, se = run_cmd(cmd, retries=args.retries, retry_sleep=args.retry_sleep,
                         timeout=getattr(args, "blast_timeout", 0), output_path=part_file)
    shutil.rmtree(tmp_dir, ignore_errors=True)
    if rc != 0:
        logger.error(f"[MMSEQS_FAIL] {batch.batch_id}\n{se[-2000:]}")
        return (batch.batch_id, False, out_file, col_names)
    part_file.touch()
    _rewrite_with_qcovs(part_file)
    os.replace(part_file, out_file)
    cache.commit(out_file, signature)
    return (batch.batch_id, True, out_file, col_names)
