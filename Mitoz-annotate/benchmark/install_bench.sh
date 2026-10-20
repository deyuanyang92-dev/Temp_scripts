#!/usr/bin/env bash
# install_bench.sh - measure how hard each annotation tool is to install from scratch.
#
# Linux: each tool is installed in a fresh ubuntu:24.04 container following its own instructions
# (knowledge/tools/*.md); wall time, conda environments, disk use and the result are recorded.
# macOS: only a solver check (mamba dry-run with CONDA_SUBDIR osx-64 / osx-arm64), marked "dry-run (未实测)".
# Appends rows to benchmark/sets/install_<YYYY-MM>.tsv:
#   tool platform method steps manual_downloads conda_envs disk_gb minutes result notes
#
# Usage: bash install_bench.sh [ours mitoz mitos2 mitofinder]   (default: all)
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPE="$(dirname "$HERE")"
OUT="$HERE/sets/install_$(date +%Y-%m).tsv"
LOGDIR="${INSTALL_BENCH_LOGS:-$PIPE/../../projects/mito_benchmark/install_logs}"
TOOLS="${*:-ours mitoz mitos2 mitofinder}"
mkdir -p "$LOGDIR"
[ -s "$OUT" ] || printf 'tool\tplatform\tmethod\tsteps\tmanual_downloads\tconda_envs\tdisk_gb\tminutes\tresult\tnotes\n' > "$OUT"

# Miniforge bootstrap shared by the per-tool recipes (counts as one step for every tool except ours,
# whose install.sh does it itself)
BOOT='apt-get update -qq && apt-get install -y -qq curl bzip2 ca-certificates >/dev/null &&
curl -sSL -o /tmp/mf.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh &&
bash /tmp/mf.sh -b -p /opt/conda >/dev/null && export PATH=/opt/conda/bin:$PATH'

# tool -> "steps|manual_downloads|recipe" (recipe runs inside the container)
recipe() {
  case "$1" in
    ours) echo '2|0|apt-get update -qq && apt-get install -y -qq curl bzip2 ca-certificates python3 >/dev/null &&
      bash /src/install.sh --refdir /opt/mitos_refs && ls ~/miniforge3/envs /opt/conda/envs 2>/dev/null' ;;
    mitoz) echo "2|0|$BOOT && mamba create -y -q -n mitoz -c conda-forge -c bioconda mitoz=3.6 >/dev/null &&
      /opt/conda/envs/mitoz/bin/mitoz --version" ;;
    mitos2) echo "3|1|$BOOT && mamba create -y -q -n mitos -c conda-forge -c bioconda mitos=2.1.10 >/dev/null &&
      mkdir -p /opt/mitos_refs && curl -sSL https://zenodo.org/records/4284483/files/refseq89m.tar.bz2 | tar -xj -C /opt/mitos_refs &&
      PATH=/opt/conda/envs/mitos/bin:\$PATH runmitos --version" ;;
    mitofinder) echo "2|0|$BOOT && mamba create -y -q -n mitofinder -c conda-forge -c bioconda mitofinder >/dev/null &&
      /opt/conda/envs/mitofinder/bin/mitofinder --version" ;;
    *) return 1 ;;
  esac
}

pkg() {   # conda package spec for the macOS solver check
  case "$1" in
    ours) echo "mitoz=3.6 exonerate blast infernal table2asn aragorn trnascan-se" ;;   # install.sh skips MITOS2 on macOS
    mitoz) echo "mitoz=3.6" ;; mitos2) echo "mitos=2.1.10" ;; mitofinder) echo "mitofinder" ;;
  esac
}

for t in $TOOLS; do
  IFS='|' read -r steps manual cmd <<<"$(recipe "$t" | tr '\n' ' ')"
  log="$LOGDIR/$t.linux.log"
  t0=$(date +%s)
  docker run --rm -v "$PIPE":/src:ro ubuntu:24.04 bash -c "set -e; { $cmd; } || exit 9; echo ENVS \$(ls -d /opt/conda/envs/* ~/miniforge3/envs/* 2>/dev/null | wc -l); echo DISK \$(du -sk /opt/conda ~/miniforge3 /opt/mitos_refs 2>/dev/null | awk '{s+=\$1} END {print s}')" > "$log" 2>&1
  rc=$?
  mins=$(awk -v s=$(( $(date +%s) - t0 )) 'BEGIN{printf "%.1f", s/60}')
  envs=$(grep -o '^ENVS [0-9]*' "$log" | awk '{print $2}')
  disk=$(grep -o '^DISK [0-9]*' "$log" | awk '{printf "%.1f", $2/1048576}')
  res=ok; [ $rc -ne 0 ] && res="failed: $(tail -1 "$log" | tr '\t' ' ' | cut -c1-120)"
  printf '%s\tlinux-x86_64\tdocker ubuntu:24.04 fresh\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$t" "$steps" "$manual" "${envs:-}" "${disk:-}" "$mins" "$res" "log $log" >> "$OUT"
  for sub in osx-64 osx-arm64; do
    dlog="$LOGDIR/$t.$sub.log"
    CONDA_SUBDIR=$sub CONDA_OVERRIDE_OSX=13.0 mamba create --dry-run -q -n probe_$$ -c conda-forge -c bioconda $(pkg "$t") > "$dlog" 2>&1
    if [ $? -eq 0 ]; then r=solvable; else r="unsolvable: $(grep -m1 -iE 'nothing provides|does not exist|conflict|not available' "$dlog" | tr '\t' ' ' | cut -c1-120)"; fi
    note="log $dlog"; [ "$sub" = osx-arm64 ] && note="$note; install.sh uses osx-64 packages under Rosetta 2 on Apple Silicon"
    printf '%s\t%s\tdry-run (未实测)\t\t\t\t\t\t%s\t%s\n' "$t" "$sub" "$r" "$note" >> "$OUT"
  done
  echo "[install_bench] $t done"
done
echo "[install_bench] -> $OUT"
