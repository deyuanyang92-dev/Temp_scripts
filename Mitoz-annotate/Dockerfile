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
ENV PATH=/root/.local/bin:/root/miniforge3/bin:$PATH
WORKDIR /data
CMD ["mitom", "doctor"]
