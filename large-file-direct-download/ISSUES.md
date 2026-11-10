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
