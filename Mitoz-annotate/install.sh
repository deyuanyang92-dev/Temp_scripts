#!/usr/bin/env bash
# install.sh - one command installs everything the mitogenome annotation workflow needs:
#   conda environments mitoz (MitoZ 3.6), mitos (MITOS2 2.1.10), mitotools (exonerate, BLAST, Infernal,
#   table2asn, ARAGORN, tRNAscan-SE, HMMER, MAFFT), the MITOS2 reference data refseq89m (Zenodo 4284483, md5 checked),
#   Biopython for the system python3, the `mitom` command (~/.local/bin/mitom), then mito_env.py doctor
#   --write-config so no tool path has to be typed.
# Existing environments are kept (--force recreates them). Safe to run again.
#
# Platforms (checked with mamba solves, 2026-10-20):
#   Linux x86_64 / Windows WSL2   everything
#   macOS Intel                   MitoZ + mitotools; MITOS2 has no macOS build (viennarna<2) -> MitoZ's own
#                                 tRNA step is used automatically (slower). Full workflow: use the Docker image.
#   macOS Apple Silicon           as macOS Intel, through Rosetta 2 (osx-64 packages)
#   Windows (native)              not supported by the bioinformatics packages: use WSL2 or Docker
#
# Usage: bash install.sh [--force] [--refdir DIR] [--envs "mitoz mitos mitotools"]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REFDIR="${MITOS_REFDIR:-$HOME/mitos_refs}"
ENVS="mitoz mitos mitotools"
FORCE=0
REF_URL="https://zenodo.org/records/4284483/files/refseq89m.tar.bz2"
REF_MD5="bb6325e27612e61a6995e88e8ffcecf8"

while [ $# -gt 0 ]; do
  case "$1" in
    --force) FORCE=1 ;;
    --refdir) REFDIR="$2"; shift ;;
    --envs) ENVS="$2"; shift ;;
    -h|--help) sed -n '2,16p' "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done

say() { printf '\n[install] %s\n' "$*"; }

# ── platform ─────────────────────────────────────────────────────────────────
OS="$(uname -s)"; ARCH="$(uname -m)"
case "$OS" in
  Linux) ;;
  Darwin)
    if [ "$ARCH" = "arm64" ]; then
      if ! arch -x86_64 /usr/bin/true 2>/dev/null; then
        echo "Apple Silicon: Rosetta 2 is needed. Run: softwareupdate --install-rosetta --agree-to-license" >&2
        exit 1
      fi
      export CONDA_SUBDIR=osx-64          # bioconda has no arm64 build of tbl2asn (MitoZ) / table2asn
      say "Apple Silicon: installing Intel (osx-64) packages, run through Rosetta 2"
    fi
    keep=""; for e in $ENVS; do [ "$e" != mitos ] && keep="$keep $e"; done; ENVS="${keep# }"
    say "macOS: MITOS2 cannot be installed (no macOS build of viennarna<2); MitoZ's own tRNA step will be used"
    ;;
  MINGW*|MSYS*|CYGWIN*)
    echo "Windows: please use WSL2 (wsl --install -d Ubuntu, then run this script inside Ubuntu) or Docker." >&2
    exit 1 ;;
  *) echo "unsupported system: $OS" >&2; exit 1 ;;
esac

# ── conda / mamba (Miniforge installed into ~/miniforge3 when none is present) ─
SOLVER="$(command -v mamba || command -v micromamba || command -v conda || true)"
if [ -z "$SOLVER" ]; then
  for d in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" /opt/conda; do
    for x in mamba conda; do [ -z "$SOLVER" ] && [ -x "$d/bin/$x" ] && SOLVER="$d/bin/$x"; done
  done
fi
if [ -z "$SOLVER" ]; then
  say "no conda/mamba found: installing Miniforge into $HOME/miniforge3"
  MF_OS="$OS"; MF_ARCH="$ARCH"; [ -n "${CONDA_SUBDIR:-}" ] && MF_ARCH="x86_64"
  curl -fsSL -o /tmp/miniforge.sh \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-${MF_OS}-${MF_ARCH}.sh"
  bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
  rm -f /tmp/miniforge.sh
  SOLVER="$HOME/miniforge3/bin/mamba"
fi
say "using $SOLVER"

env_exists() { "$SOLVER" env list 2>/dev/null | awk '{print $1}' | grep -qx "$1"; }

for e in $ENVS; do
  if env_exists "$e" && [ "$FORCE" = 0 ]; then
    say "environment '$e' exists - kept (use --force to recreate)"
    continue
  fi
  say "creating environment '$e' from env/$e.yml"
  if env_exists "$e"; then "$SOLVER" env remove -y -n "$e"; fi
  "$SOLVER" env create -y -n "$e" -f "$HERE/env/$e.yml" 2>/dev/null || "$SOLVER" env create -n "$e" -f "$HERE/env/$e.yml"
done

# ── MITOS2 reference data ────────────────────────────────────────────────────
if [[ " $ENVS " == *" mitos "* ]] || env_exists mitos; then
  if [ -d "$REFDIR/refseq89m" ]; then
    say "MITOS2 reference data present: $REFDIR/refseq89m"
  else
    say "downloading MITOS2 reference data refseq89m (19.6 MB) to $REFDIR"
    mkdir -p "$REFDIR"
    curl -fL --retry 3 -o "$REFDIR/refseq89m.tar.bz2" "$REF_URL"
    sum="$( (md5sum "$REFDIR/refseq89m.tar.bz2" 2>/dev/null || md5 -r "$REFDIR/refseq89m.tar.bz2") | awk '{print $1}')"
    if [ "$sum" != "$REF_MD5" ]; then
      echo "checksum mismatch for refseq89m.tar.bz2 ($sum, expected $REF_MD5) - download again" >&2
      rm -f "$REFDIR/refseq89m.tar.bz2"; exit 1
    fi
    tar -xjf "$REFDIR/refseq89m.tar.bz2" -C "$REFDIR"
  fi
fi

# ── Biopython for the system python3 (runs batch_mitoz.py and the other scripts) ─
if ! python3 -c "import Bio" 2>/dev/null; then
  say "installing Biopython for $(command -v python3)"
  python3 -m pip install --user biopython 2>/dev/null || python3 -m pip install --user --break-system-packages biopython \
    || echo "could not install Biopython; run the scripts with: $SOLVER run -n mitotools python3 ..." >&2
fi

# ── the mitom command ────────────────────────────────────────────────────────
mkdir -p "$HOME/.local/bin"
ln -sf "$HERE/mitom" "$HOME/.local/bin/mitom"
chmod +x "$HERE/mitom"
case ":$PATH:" in
  *":$HOME/.local/bin:"*) say "command installed: mitom" ;;
  *) say "command installed: $HOME/.local/bin/mitom - add it to PATH: echo 'export PATH=\$HOME/.local/bin:\$PATH' >> ~/.bashrc   (zsh, the macOS default shell: ~/.zshrc instead of ~/.bashrc)" ;;
esac

# ── record what was found ────────────────────────────────────────────────────
say "checking"
MITOS_REFDIR="$REFDIR" python3 "$HERE/mito_env.py" doctor --write-config
