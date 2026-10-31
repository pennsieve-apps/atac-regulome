# --- MACS3 has no wheels on PyPI: build one (pinned compiler and build tools) -------------------
FROM --platform=linux/amd64 python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337 AS macs3
RUN apt-get update \
    && apt-get install -y --no-install-recommends gcc=4:12.2.0-3 libc6-dev=2.36-9+deb12u14 zlib1g-dev=1:1.2.13.dfsg-1 \
    && rm -rf /var/lib/apt/lists/*
ADD --checksum=sha256:22cb5e22903ef0598e8078cc19d4593cb43be43870808335c98e0f7298597224 \
    https://files.pythonhosted.org/packages/b2/4b/d228c770577e48a11bdad778130590f0f68e7659e3d5bafa1bdef7ce25ad/macs3-3.0.5.tar.gz /src/
RUN pip install --no-cache-dir setuptools==84.2.1 numpy==2.5.4 Cython==3.3.0 \
    && pip wheel --no-cache-dir --no-deps --no-build-isolation -w /wheels /src/macs3-3.0.5.tar.gz

# --- Reference data: GENCODE annotations + genomes (the files SnapATAC2 2.10.0's genome objects
# download lazily, with the sha256 SnapATAC2 pins for them) and JASPAR 2026 CORE vertebrates -----
FROM --platform=linux/amd64 python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337 AS reference
RUN pip install --no-cache-dir pyfaidx==0.8.2 numpy==2.5.4 snapatac2==2.10.0 --no-deps \
    && pip install --no-cache-dir anndata==0.13.4 pandas==2.3.3 polars==1.44.2 polars-runtime-32==1.44.2 scipy==1.18.1 \
       pooch==1.9.0 typeguard==4.6.0 scverse-misc==0.1.8 rustworkx==0.18.1 scikit-learn==1.9.1 natsort==8.4.0 h5py==3.16.0 \
       pyarrow==26.0.0 tqdm==4.70.1
ADD --checksum=sha256:b82a655bdb736ca0e463a8f5d00242bedf10fa88ce9d651a017f135c7c4e9285 \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.basic.annotation.gff3.gz \
    /opt/genomes/GRCh38/annotation.gff3.gz
ADD --checksum=sha256:6f433e2676e26569a678ce78b37e94a64ddd50a09479e433ad6f75e37dc82e48 \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/gencode.vM30.basic.annotation.gff3.gz \
    /opt/genomes/GRCm39/annotation.gff3.gz
ADD --checksum=sha256:4fac949d7021cbe11117ddab8ec1960004df423d672446cadfbc8cca8007e228 \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz /src/GRCh38.fa.gz
ADD --checksum=sha256:3b923c06a0d291fe646af6bf7beaed7492bf0f6dd5309d4f5904623cab41b0aa \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/GRCm39.primary_assembly.genome.fa.gz /src/GRCm39.fa.gz
ADD --checksum=sha256:ae4f77a17256c7b8b7b73402810c1ecf4fc411e9ad7832cdc87d586dee9a4a15 \
    https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt /opt/motifs/
COPY scripts/primary_fasta.py /src/
RUN for g in GRCh38 GRCm39; do \
        python -c "import snapatac2 as s, sys; sys.stdout.write(''.join(f'{c}\t{n}\n' for c, n in getattr(s.genome, '$g').chrom_sizes.items()))" > /opt/genomes/$g/chrom_sizes.tsv \
        && python /src/primary_fasta.py /src/$g.fa.gz /opt/genomes/$g/chrom_sizes.tsv /opt/genomes/$g \
        && rm /src/$g.fa.gz; \
    done && chmod -R a+rX /opt/genomes /opt/motifs

# --- The app -------------------------------------------------------------------------------------
FROM --platform=linux/amd64 python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /processor

# Pinned Python dependencies, binary wheels only, installed without resolving (requirements.txt is
# the complete set; see its header), plus the MACS3 wheel built above.
COPY requirements.txt /processor/requirements.txt
COPY --from=macs3 /wheels /wheels
RUN pip install --only-binary=:all: --no-deps -r /processor/requirements.txt /wheels/macs3-3.0.5-*.whl \
    && rm -rf /wheels \
    && find /usr/local/lib/python3.12/site-packages -name '*.pyc' -delete \
    && find /usr/local/lib/python3.12/site-packages -depth -type d \( -name tests -o -name test \) -exec rm -rf {} + \
    && rm -rf /usr/local/lib/python3.12/site-packages/pyarrow/include

COPY --from=reference /opt/genomes /opt/genomes
COPY --from=reference /opt/motifs /opt/motifs

COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH="/processor" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    GENOMES_DIR=/opt/genomes \
    MOTIFS_FILE=/opt/motifs/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib \
    NUMBA_CACHE_DIR=/tmp/numba \
    HOME=/tmp

# Imports checked and matplotlib's font cache built now; the cache folders left in /tmp are made writable
# for uid 1000 (the platform does not remount /tmp).
RUN python -c "import matplotlib.pyplot, snapatac2, MACS3.Signal.PeakDetect, igraph" \
    && mkdir -p /tmp/matplotlib /tmp/numba && chmod -R a+rwX /tmp/matplotlib /tmp/numba

# Runs as uid 1000 on the platform; nothing here needs root. Writes only to OUTPUT_DIR and /tmp.
USER 1000

ENTRYPOINT ["python3", "-m", "processor.entry"]
