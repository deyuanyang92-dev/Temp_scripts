"""`blast2metadata-db`: find, inspect, update and clean NCBI BLAST databases (nt, nr, ...).

NCBI publishes one JSON manifest per database (https://ftp.ncbi.nlm.nih.gov/blast/db/<db>-<nucl|prot>-metadata.json)
listing the tarballs, the volume count and `last-updated`; each tarball has a `.md5` next to it.
`update` compares that manifest with the local database date (blastdbcmd -info), downloads the tarballs
(aria2c with md5 checksums and retry rounds, or ascp, or NCBI's update_blastdb.pl), verifies them,
extracts into a staging directory, checks the new database and then swaps it in.

The downloaded BLAST v5 databases are used directly by BLAST+, by DIAMOND >= 2.1.14 (protein DBs)
and by `blast2metadata --engine blast|diamond`.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tarfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.request import Request, urlopen

from . import __version__

NCBI_DB_URL = "https://ftp.ncbi.nlm.nih.gov/blast/db/"
ALIAS_SUFFIX = {"nucl": ".nal", "prot": ".pal"}
INDEX_SUFFIX = {"nucl": ".nin", "prot": ".pin"}
VOLUME_RE = re.compile(r"^(?P<name>.+?)\.(?P<vol>\d{2,3})\.(?P<ext>[np](?:in|hr|sq|og|os|ot|tf|to|db|al))$")


def _log(msg: str) -> None:
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def http_get(url: str, timeout: int = 60) -> bytes:
    req = Request(url, headers={"User-Agent": f"blast2metadata-db/{__version__}"})
    with urlopen(req, timeout=timeout) as r:
        return r.read()


# ─────────────────────────────────────────────
# local databases
# ─────────────────────────────────────────────
def blastdb_dirs(extra: List[str]) -> List[Path]:
    dirs: List[Path] = []
    for d in list(extra) + [x for x in os.environ.get("BLASTDB", "").split(os.pathsep) if x]:
        p = Path(d).expanduser()
        if p.is_dir() and p.resolve() not in [x.resolve() for x in dirs]:
            dirs.append(p)
    return dirs


def find_local_dbs(dirs: List[Path]) -> List[Tuple[Path, str]]:
    """(base_path, 'nucl'|'prot') for every database in *dirs* (alias file, or single-volume index)."""
    found: List[Tuple[Path, str]] = []
    for d in dirs:
        names = set()
        for f in d.iterdir():
            for t, suf in ALIAS_SUFFIX.items():
                if f.name.endswith(suf) and not VOLUME_RE.match(f.name):
                    names.add((f.name[: -len(suf)], t))
            for t, suf in INDEX_SUFFIX.items():
                if f.name.endswith(suf) and not VOLUME_RE.match(f.name):
                    names.add((f.name[: -len(suf)], t))
        for name, t in sorted(names):
            found.append((d / name, t))
    return found


def blastdbcmd_info(base: Path, dbtype: str) -> Dict[str, object]:
    """Parse `blastdbcmd -info`: sequences, letters, date (as 'YYYY-MM-DD'), version, volumes."""
    exe = shutil.which("blastdbcmd")
    info: Dict[str, object] = {"ok": False}
    if not exe:
        info["error"] = "blastdbcmd not in PATH"
        return info
    p = subprocess.run([exe, "-db", str(base), "-dbtype", dbtype, "-info"], capture_output=True, text=True)
    if p.returncode != 0:
        info["error"] = (p.stderr or p.stdout).strip()[-300:]
        return info
    txt = p.stdout
    m = re.search(r"([\d,]+) sequences; ([\d,]+) total (?:bases|residues)", txt)
    if m:
        info["sequences"] = int(m.group(1).replace(",", ""))
        info["letters"] = int(m.group(2).replace(",", ""))
    m = re.search(r"Date: (\w{3}) (\d{1,2}), (\d{4})", txt)
    if m:
        months = {m_: i for i, m_ in enumerate("Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split(), 1)}
        info["date"] = f"{m.group(3)}-{months[m.group(1)]:02d}-{int(m.group(2)):02d}"
    m = re.search(r"BLASTDB Version: (\d+)", txt)
    if m:
        info["version"] = int(m.group(1))
    info["volumes"] = len(re.findall(r"^\s+\S+\.\d{3}$", txt, flags=re.M)) or 1
    info["ok"] = True
    return info


def leftovers(d: Path, name: str) -> Tuple[List[Path], int]:
    files = sorted(list(d.glob(f"{name}.*.tar.gz")) + list(d.glob(f"{name}.tar.gz")) + list(d.glob(f"{name}.*.aria2"))
                   + list(d.glob(f"{name}.*.tar.gz.md5")) + list(d.glob(f"{name}.tar.gz.md5")))
    return files, sum(f.stat().st_size for f in files if f.is_file())


# ─────────────────────────────────────────────
# remote manifest
# ─────────────────────────────────────────────
def fetch_manifest(name: str, dbtype: Optional[str] = None) -> Dict[str, object]:
    types = [dbtype] if dbtype else ["nucl", "prot"]
    last: Optional[Exception] = None
    for t in types:
        try:
            d = json.loads(http_get(f"{NCBI_DB_URL}{name}-{t}-metadata.json").decode("utf-8"))
            d["_dbtype"] = t
            d["files"] = [u.replace("ftp://", "https://") for u in d.get("files", [])]
            d["_last_updated"] = str(d.get("last-updated", ""))[:10]
            return d
        except Exception as e:  # 404 for the wrong type
            last = e
    raise SystemExit(f"ERROR: no NCBI manifest for database '{name}' ({NCBI_DB_URL}{name}-nucl|prot-metadata.json): {last}")


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{int(n)} B"
        n /= 1024
    return f"{n:.1f} TB"


# ─────────────────────────────────────────────
# download
# ─────────────────────────────────────────────
def fetch_md5s(urls: List[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for u in urls:
        try:
            txt = http_get(u + ".md5").decode("utf-8", "replace").split()
            if txt:
                out[Path(u).name] = txt[0].lower()
        except Exception as e:
            _log(f"[MD5] could not fetch {u}.md5: {e}")
    return out


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(16 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_aria2(urls: List[str], md5s: Dict[str, str], outdir: Path, jobs: int, split: int, conn: int,
                   rounds: int, aria2: str = "aria2c") -> List[str]:
    """Download with aria2c, letting it verify md5 checksums; retry the failed ones for *rounds*.
    Returns the file names that are still missing/bad."""
    exe = shutil.which(aria2)
    if not exe:
        raise SystemExit(f"ERROR: aria2c not found ({aria2}); use --downloader update_blastdb or install aria2")
    outdir.mkdir(parents=True, exist_ok=True)
    pending = list(urls)
    for rnd in range(1, rounds + 1):
        todo = [u for u in pending if not _tarball_ok(outdir / Path(u).name, md5s.get(Path(u).name))]
        if not todo:
            return []
        lst = outdir / f".aria2_round{rnd}.list"
        with open(lst, "w", encoding="utf-8") as fh:
            for u in todo:
                fh.write(u + "\n")
                m = md5s.get(Path(u).name)
                if m:
                    fh.write(f"  checksum=md5={m}\n")
        cmd = [exe, "--continue=true", "--allow-overwrite=true", "--auto-file-renaming=false",
               "--file-allocation=trunc", f"--max-concurrent-downloads={jobs}", f"--split={split}",
               f"--max-connection-per-server={conn}", "--min-split-size=50M", "--max-tries=10",
               "--retry-wait=10", "--connect-timeout=30", "--timeout=60", "--summary-interval=30",
               "--console-log-level=warn", "-d", str(outdir), "-i", str(lst)]
        _log(f"[ARIA2] round {rnd}/{rounds}: {len(todo)} file(s)  ->  {' '.join(cmd[:1])} ... -i {lst.name}")
        # aria2 only speaks HTTP proxies: a socks:// all_proxy (common with Clash/V2Ray setups) makes it
        # abort at start-up, so drop it and keep http_proxy/https_proxy.
        env = {k: v for k, v in os.environ.items() if k.lower() != "all_proxy"}
        rc = subprocess.run(cmd, env=env).returncode
        pending = todo
        if rc != 0 and exe.startswith("/snap/") and not any((outdir / Path(u).name).exists() for u in todo):
            _log("[ARIA2] hint: the snap-packaged aria2c cannot read /tmp; use a destination under $HOME or /media, "
                 "or install aria2 with apt/conda")
        if rc != 0 and rnd < rounds:
            time.sleep(10)
    return [Path(u).name for u in pending if not _tarball_ok(outdir / Path(u).name, md5s.get(Path(u).name))]


def _tarball_ok(path: Path, md5: Optional[str]) -> bool:
    if not path.is_file() or path.stat().st_size == 0 or path.with_name(path.name + ".aria2").exists():
        return False
    if md5:
        ok = md5_of(path) == md5
        if not ok:
            _log(f"[MD5] mismatch, will re-download: {path.name}")
        return ok
    return True


def download_ascp(urls: List[str], md5s: Dict[str, str], outdir: Path, jobs: int, rate: str, key: str) -> List[str]:
    exe = shutil.which("ascp")
    if not exe:
        raise SystemExit("ERROR: ascp (Aspera Connect) not found in PATH")
    key = os.path.expanduser(key)
    if not os.path.isfile(key):
        raise SystemExit(f"ERROR: Aspera key not found: {key}")
    outdir.mkdir(parents=True, exist_ok=True)
    import concurrent.futures as cf

    def one(u: str) -> bool:
        name = Path(u).name
        if _tarball_ok(outdir / name, md5s.get(name)):
            return True
        cmd = [exe, "-T", "-k", "1", "-Q", "-l", rate, "-i", key,
               f"anonftp@ftp-private.ncbi.nlm.nih.gov:/blast/db/{name}", str(outdir)]
        for attempt in range(1, 4):
            _log(f"[ASCP] {name} attempt {attempt}")
            if subprocess.run(cmd, capture_output=True).returncode == 0 and _tarball_ok(outdir / name, md5s.get(name)):
                return True
        return False

    with cf.ThreadPoolExecutor(max_workers=max(1, jobs)) as ex:
        results = list(ex.map(one, urls))
    return [Path(u).name for u, ok in zip(urls, results) if not ok]


def download_update_blastdb(name: str, outdir: Path, exe: str) -> None:
    """Let NCBI's update_blastdb.pl do download + md5 check + decompress inside *outdir*."""
    script = shutil.which(exe) or exe
    if not os.path.isfile(script):
        raise SystemExit(f"ERROR: update_blastdb.pl not found: {exe}")
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = ["perl", script, "--decompress", "--blastdb_version", "5", "--num_threads", "4", name]
    _log(f"[UPDATE_BLASTDB] {' '.join(cmd)}  (cwd={outdir})")
    p = subprocess.run(cmd, cwd=str(outdir))
    if p.returncode != 0:
        raise SystemExit(f"ERROR: update_blastdb.pl failed (rc={p.returncode})")


# ─────────────────────────────────────────────
# extract + swap
# ─────────────────────────────────────────────
def extract_all(tarballs: List[Path], staging: Path) -> None:
    staging.mkdir(parents=True, exist_ok=True)
    for i, tb in enumerate(tarballs, 1):
        _log(f"[EXTRACT] {i}/{len(tarballs)} {tb.name}")
        with tarfile.open(tb, "r:gz") as tf:
            for m in tf.getmembers():  # NCBI tarballs are flat; refuse anything that escapes the directory
                if m.name.startswith("/") or ".." in Path(m.name).parts:
                    raise SystemExit(f"ERROR: unsafe path in {tb.name}: {m.name}")
            tf.extractall(staging)


def swap_in(staging: Path, dest: Path, name: str, dbtype: str, keep_old: bool) -> Path:
    """Move the new database files from *staging* into *dest*; old files of the same database go to
    dest/<name>_old_<date>/ (or are deleted)."""
    old_dir = dest / f"{name}_old_{time.strftime('%Y%m%d')}"
    old_files = [f for f in dest.iterdir() if f.is_file() and (
        f.name == name + ALIAS_SUFFIX[dbtype] or f.name.startswith(name + ".") and not f.name.endswith((".tar.gz", ".md5", ".aria2")))]
    if old_files:
        if keep_old:
            old_dir.mkdir(exist_ok=True)
            for f in old_files:
                shutil.move(str(f), str(old_dir / f.name))
            _log(f"[SWAP] moved {len(old_files)} old file(s) to {old_dir}")
        else:
            for f in old_files:
                f.unlink()
            _log(f"[SWAP] deleted {len(old_files)} old file(s)")
    n = 0
    for f in sorted(staging.iterdir()):
        if f.name.endswith((".md5", ".aria2")) or f.name.startswith(".aria2_round"):
            f.unlink()
            continue
        shutil.move(str(f), str(dest / f.name))
        n += 1
    staging.rmdir()
    _log(f"[SWAP] installed {n} file(s) into {dest}")
    return dest / name


# ─────────────────────────────────────────────
# commands
# ─────────────────────────────────────────────
def cmd_list(args: argparse.Namespace) -> int:
    dirs = blastdb_dirs(args.dir)
    if not dirs:
        print("No BLAST database directory found: set $BLASTDB or pass --dir")
        return 1
    dbs = find_local_dbs(dirs)
    if not dbs:
        print(f"No BLAST databases in: {', '.join(map(str, dirs))}")
        return 1
    for base, t in dbs:
        info = blastdbcmd_info(base, t)
        line = f"{base}  [{t}]"
        if info.get("ok"):
            line += (f"  v{info.get('version', '?')}  {info.get('date', '?')}  {info.get('sequences', 0):,} seqs  "
                     f"{info.get('volumes', 1)} vol")
        else:
            line += f"  ({info.get('error', 'unreadable')})"
        files, size = leftovers(base.parent, base.name)
        if files:
            line += f"  leftovers: {len(files)} tarball/.aria2 files = {human(size)} (clean with `clean`)"
        print(line)
        if args.remote and info.get("ok"):
            try:
                man = fetch_manifest(base.name, t)
            except SystemExit as e:
                print(f"    remote: {e}")
                continue
            newer = man["_last_updated"] > str(info.get("date", ""))
            print(f"    remote: {man['_last_updated']}  {man.get('number-of-sequences', 0):,} seqs  "
                  f"{man.get('number-of-volumes', 0)} vol  download {human(man.get('bytes-total-compressed', 0))}  "
                  f"on disk {human(man.get('bytes-total', 0))}  ->  "
                  + ("UPDATE AVAILABLE" if newer else "up to date"))
    return 0


def cmd_update(args: argparse.Namespace) -> int:
    name = args.name
    dest = Path(args.dir).expanduser().resolve()
    dest.mkdir(parents=True, exist_ok=True)
    man = fetch_manifest(name, args.dbtype)
    dbtype = man["_dbtype"]
    urls: List[str] = list(man["files"])
    _log(f"[REMOTE] {name} ({dbtype}) last-updated {man['_last_updated']}: {len(urls)} tarball(s), "
         f"{human(man.get('bytes-total-compressed', 0))} to download, {human(man.get('bytes-total', 0))} on disk")
    local = blastdbcmd_info(dest / name, dbtype) if (dest / name).with_name(name + ALIAS_SUFFIX[dbtype]).exists() or \
        (dest / name).with_name(name + INDEX_SUFFIX[dbtype]).exists() else {"ok": False}
    if local.get("ok"):
        _log(f"[LOCAL] {dest / name}: {local.get('date')} {local.get('sequences', 0):,} seqs {local.get('volumes')} vol")
        if not args.force and str(local.get("date", "")) >= man["_last_updated"]:
            _log("[LOCAL] already up to date (use --force to re-download)")
            return 0
    else:
        _log(f"[LOCAL] no {name} database in {dest} yet")
    free = shutil.disk_usage(dest).free
    need = int(man.get("bytes-total-compressed", 0)) + int(man.get("bytes-total", 0))
    if free < need:
        _log(f"[DISK] WARNING: {human(free)} free in {dest}, download + extraction need about {human(need)}")
        if not args.force:
            raise SystemExit("ERROR: not enough disk space (use --force to try anyway, or --keep_tarballs off and --in_place)")
    if args.dry_run:
        print("\n".join(urls))
        return 0

    dl_dir = dest / f".download_{name}"
    if args.downloader == "update_blastdb":
        download_update_blastdb(name, dl_dir, args.update_blastdb_exe)
        staging = dl_dir  # it already decompressed in place
        tarballs = []
    else:
        md5s = fetch_md5s(urls)
        if args.downloader == "aria2":
            bad = download_aria2(urls, md5s, dl_dir, args.jobs, args.split, args.conn, args.rounds, args.aria2_exe)
        else:
            bad = download_ascp(urls, md5s, dl_dir, args.jobs, args.ascp_rate, args.ascp_key)
        if bad:
            raise SystemExit(f"ERROR: {len(bad)} tarball(s) still missing or failed md5 after retries: {bad[:10]} "
                             f"(rerun the same command to resume)")
        tarballs = [dl_dir / Path(u).name for u in urls]
        staging = dest / f".staging_{name}"
        if staging.exists():
            shutil.rmtree(staging)
        extract_all(tarballs, staging)
    new_info = blastdbcmd_info(staging / name, dbtype)
    if not new_info.get("ok"):
        raise SystemExit(f"ERROR: extracted database fails blastdbcmd -info: {new_info.get('error')} (left in {staging})")
    _log(f"[NEW] {new_info.get('date')} {new_info.get('sequences', 0):,} seqs {new_info.get('volumes')} vol  v{new_info.get('version')}")
    base = swap_in(staging, dest, name, dbtype, keep_old=args.keep_old)
    if not args.keep_tarballs:
        for tb in tarballs:
            tb.unlink(missing_ok=True)
        for f in dl_dir.glob("*"):
            if f.suffix in {".list", ".aria2"} or f.name.startswith(".aria2_round"):
                f.unlink(missing_ok=True)
        if dl_dir.exists() and not any(dl_dir.iterdir()):
            dl_dir.rmdir()
    with open(dest / f"{name}.update.json", "w", encoding="utf-8") as fh:
        json.dump({"name": name, "dbtype": dbtype, "remote": {k: v for k, v in man.items() if k != "files"},
                   "installed": new_info, "when": time.strftime("%Y-%m-%d %H:%M:%S"), "downloader": args.downloader},
                  fh, indent=2, default=str)
    _log(f"[DONE] {base}  (export BLASTDB={dest})")
    return 0


def cmd_clean(args: argparse.Namespace) -> int:
    d = Path(args.dir).expanduser().resolve()
    files, size = leftovers(d, args.name)
    for sub in d.glob(f".download_{args.name}"):
        files += sorted(sub.glob("*"))
    if not files:
        print(f"nothing to clean for {args.name} in {d}")
        return 0
    print(f"{len(files)} leftover file(s), {human(size)}:")
    for f in files[:15]:
        print("  ", f.name)
    if len(files) > 15:
        print(f"   ... {len(files) - 15} more")
    if not args.yes:
        print("re-run with --yes to delete them")
        return 0
    for f in files:
        f.unlink(missing_ok=True)
    print("deleted")
    return 0


def cmd_env(args: argparse.Namespace) -> int:
    dirs = blastdb_dirs(args.dir)
    if not dirs:
        print("No BLAST database directory found: pass --dir")
        return 1
    print("export BLASTDB=" + os.pathsep.join(str(d.resolve()) for d in dirs))
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="blast2metadata-db",
                                description="Find, inspect, update and clean NCBI BLAST databases (nt, nr, swissprot, ...).")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sp = p.add_subparsers(dest="cmd", required=True)

    a = sp.add_parser("list", help="databases in $BLASTDB / --dir, with --remote: compare with NCBI")
    a.add_argument("--dir", action="append", default=[], help="directory to scan (repeatable; $BLASTDB is always scanned)")
    a.add_argument("--remote", action="store_true", help="fetch the NCBI manifest and say whether an update exists")
    a.set_defaults(func=cmd_list)

    u = sp.add_parser("update", help="download/update one NCBI database into --dir (resumable)")
    u.add_argument("name", help="NCBI database name: nt, nr, core_nt, swissprot, 16S_ribosomal_RNA, ...")
    u.add_argument("-o", "--dir", required=True, help="destination directory (the database lives directly in it)")
    u.add_argument("--dbtype", choices=["nucl", "prot"], default=None, help="override manifest type lookup")
    u.add_argument("--downloader", choices=["aria2", "ascp", "update_blastdb"], default="aria2",
                   help="aria2 (default; md5-verified, resumable, honours http_proxy), ascp (Aspera), "
                        "or NCBI's update_blastdb.pl (does download + decompress itself)")
    u.add_argument("--jobs", type=int, default=4, help="parallel downloads (default 4)")
    u.add_argument("--split", type=int, default=8, help="aria2 --split per file (default 8)")
    u.add_argument("--conn", type=int, default=8, help="aria2 --max-connection-per-server (default 8)")
    u.add_argument("--rounds", type=int, default=10, help="aria2 retry rounds for failed files (default 10)")
    u.add_argument("--aria2_exe", default="aria2c")
    u.add_argument("--ascp_rate", default="200m", help="ascp -l speed limit (default 200m)")
    u.add_argument("--ascp_key", default="~/.aspera/connect/etc/asperaweb_id_dsa.openssh")
    u.add_argument("--update_blastdb_exe", default="update_blastdb.pl")
    u.add_argument("--keep_tarballs", action="store_true", help="keep the downloaded tar.gz files")
    u.add_argument("--keep_old", action="store_true", help="move the previous database to <name>_old_<date>/ instead of deleting")
    u.add_argument("--force", action="store_true", help="re-download even if the local database is up to date / disk looks short")
    u.add_argument("--dry_run", action="store_true", help="only print the tarball URLs")
    u.set_defaults(func=cmd_update)

    c = sp.add_parser("clean", help="delete leftover tar.gz / .aria2 / .md5 files of a database")
    c.add_argument("name")
    c.add_argument("-o", "--dir", required=True)
    c.add_argument("--yes", action="store_true")
    c.set_defaults(func=cmd_clean)

    e = sp.add_parser("env", help="print an `export BLASTDB=...` line for the directories found")
    e.add_argument("--dir", action="append", default=[])
    e.set_defaults(func=cmd_env)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
