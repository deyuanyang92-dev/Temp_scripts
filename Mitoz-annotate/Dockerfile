# Mitogenome annotation workflow (MitoZ + MITOS2 + mitotools + refseq89m) for any system with Docker
# (Windows, macOS incl. Apple Silicon via emulation, Linux). Built by install.sh in a clean Ubuntu.
#   docker build -t yzz-mito .
#   docker run --rm -v "$PWD":/data yzz-mito python3 /opt/mito/batch_mitoz.py run -i /data/gb -o /data/out \
#       --clade Annelida-segmented-worms --genetic_code 5
FROM --platform=linux/amd64 ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive LANG=C.UTF-8
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl bzip2 procps \
        python3 python3-biopython && rm -rf /var/lib/apt/lists/*
COPY . /opt/mito
RUN bash /opt/mito/install.sh --refdir /opt/mitos_refs && /root/miniforge3/bin/mamba clean -a -y
ENV PATH=/root/miniforge3/bin:$PATH
WORKDIR /data
CMD ["python3", "/opt/mito/mito_env.py", "doctor"]
