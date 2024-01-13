FROM mambaorg/micromamba

WORKDIR /workflow

# Install deps for dna-brnn as root
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    libz-dev \
    gcc \
    git \
    make \
    && rm -rf /var/lib/apt/lists/*

# Install dna-brnn
RUN git clone https://github.com/lh3/dna-nn && \
    cd dna-nn && \
    make && \
    # Then remove unneeded tools
    apt-get purge gcc git make -y

# Update path for dna-brnn
ENV PATH="${PATH}:/workflow/dna-nn/"

# Return to mamba user.
USER $MAMBA_USER

# Copy env.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml

# Install conda packages
RUN micromamba install --yes --file /tmp/env.yaml && \
    micromamba clean --all --yes

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "snakemake" ]
