# ISSUES / TROUBLESHOOTING

Record every run error here: date, symptom, cause, fix.

## 2026-11-10 — initial validation (FlClash, mode global, TUN off, default route eno1)

- `selftest`: curl / wget / aria2c / python each made 1 DIRECT socket to ftp.ncbi.nlm.nih.gov via eno1, 0 proxy
  sockets; negative control (curl with normal proxy env) detected as PROXY. PASSED.
- `fetch` of 2 RefSeq files (E. coli fna.gz with md5, gff.gz without): VERIFIED_CHECKSUM / VERIFIED_SIZE+gzip -t;
  re-run skipped both. Monitor: DIRECT=2, PROXY=0, egress eno1.
- Leak test 1: `run -- bash -c "https_proxy=... curl ..."` still went DIRECT — expected, child has NO_PROXY='*'
  which curl honours over https_proxy (layer-1 protection).
- Leak test 2: `run -- curl -x http://127.0.0.1:7890 --noproxy '' ...` → LEAK detected in first sample, job killed,
  exit 3 (layer-2 protection).
- Not yet tested: a real organise-sra / prefetch SRA download under `run` (inherits the clean env by design,
  organise_sra spawns kingfisher via subprocess with os.environ copy — src/organise_sra/download.py:155).

## 2026-11-10 — security review: aria2 input-file injection / path traversal in `fetch`

- Cause: `out` (from manifest, or derived from the URL with `unquote`) was written unchecked into the aria2c
  input file and joined to outdir. `%0A` in a URL or a manifest `out` like `../x` could inject aria2 options
  (e.g. `all-proxy=` → proxy use) or write outside outdir.
- Fix: every entry is validated before any network access: URL without control chars/whitespace and scheme
  http/https/ftp/sftp; `out` must be a plain file name (no `/`, `..`, control chars); md5/sha256 must be hex.
  Tested: `../evil`, `x%0A all-proxy=...`, and a manifest line `  dir=/tmp` are all refused (exit 1).
