# Mitogenome annotation workflow (MitoZ + MITOS2 + mitotools + refseq89m) for any system with Docker
# (Windows, macOS incl. Apple Silicon via emulation, Linux). Built by install.sh in a clean Ubuntu.
#   docker build -t mitom .
#   docker run --rm -v "$PWD":/data mitom mitom run -i /data/genomes.gb -o /data/out
FROM --platform=linux/amd64 ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive LANG=C.UTF-8
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl bzip2 procps \
        python3 python3-biopython && rm -rf /var/lib/apt/lists/*
COPY . /opt/mito
RUN bash /opt/mito/install.sh --refdir /opt/mitos_refs && /root/miniforge3/bin/mamba clean -a -y
# conda LAST: the system python3 (with the apt Biopython) must stay ahead of Miniforge's python3, which has none
ENV PATH=/root/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/root/miniforge3/bin
# build-time check: the launcher starts and the python that runs the scripts imports Biopython
RUN mitom run -h >/dev/null && python3 -c "import Bio"
WORKDIR /data
CMD ["mitom", "doctor"]
