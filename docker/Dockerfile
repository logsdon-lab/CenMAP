FROM condaforge/mambaforge

WORKDIR /workflow

# Install deps for dna-brnn as root
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

# Copy env.yaml
COPY env.yaml /tmp/env.yaml

# Install conda packages
ARG ENV_NAME="hgsvc3"
RUN mamba env create -f /tmp/env.yaml -n $ENV_NAME && \
    mamba clean --all --yes

# Activate environment automatically by modifying ~/.bashrc
RUN echo "conda activate $ENV_NAME" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

ENTRYPOINT [ "snakemake" ]
