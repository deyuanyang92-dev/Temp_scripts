---
name: large-file-direct-download
description: Download large research data or any GB-scale file over the local network (DIRECT) instead of the Clash/Mihomo proxy, so proxy-node traffic is not consumed. Use whenever the user asks to download SRA/ENA/GenBank/RefSeq/UniProt data, FASTQ/BAM/CRAM/FASTA, reference genomes, metagenomes, protein/BLAST/Kraken databases, GitHub release binaries over ~100 MB, or says a file is large/GB; also when running organise-sra download, prefetch, fasterq-dump, kingfisher, datasets download, or download scripts. Do not use for Claude Code's own requests, small API calls, git push/pull, or normal web fetches.
---

# Large file direct download

Tool: `direct-download` (~/.local/bin → Temp_scripts/large-file-direct-download/, README there).
It strips proxy settings only in the child process, refuses if Mihomo TUN / policy routing is active, and
proves DIRECT by watching every TCP socket of the job's process tree (kills the job on a proxy connection).

## Trigger
- Known-large data (SRA, FASTQ, BAM/CRAM, genomes, databases, archives) or the user says GB-scale → use it.
- Unknown size → `direct-download size URL` first (direct HEAD). ≥ ~100 MB → use it. Unknown → use it
  (conservative). Do not judge by extension alone.
- Never change Claude Code's own environment; never `unset` proxies in the main shell; never switch Mihomo mode,
  TUN or system routes yourself (ask the user first).

## Steps
1. `direct-download check` (once per session). Exit 2 → stop and report the PROBLEM lines; do not download via proxy.
2. Plain URLs → `direct-download fetch -o <dir> URL... [--md5 X]` or `-i manifest.tsv` (url, out, md5|sha256, size).
   Add checksums whenever the source publishes them (NCBI `md5checksums.txt`, ENA `fastq_md5`, UniProt `RELEASE.metalink`).
3. Any other tool → `direct-download run --log-dir <dir> -- <command ...>` (organise-sra, prefetch, kingfisher, datasets,
   wget, python scripts). Long jobs: run in background and read the log.
4. Read the result honestly: exit 0 only. 3 = proxy leak / DIRECT not proven, 4 = unverified integrity,
   5 = host unreachable directly (do NOT fall back to the proxy without asking), 1 = failed (rerun resumes).
   Quote the `monitor:` line (DIRECT sockets, egress dev, PROXY=0) and the report path to the user.

## Limits to tell the user when relevant
TCP only (Aspera FASP/UDP and QUIC not proven); very short connections may be missed by the 0.5 s sampler;
"DIRECT" = not via Mihomo proxy nodes, leaving via the host NIC — VPNs upstream of this host are invisible.
