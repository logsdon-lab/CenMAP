# HGSVC3
[![CI](https://github.com/logsdon-lab/hgsvc3/actions/workflows/main.yml/badge.svg)](https://github.com/logsdon-lab/hgsvc3/actions/workflows/main.yml)
[![Docker Image Version](https://img.shields.io/docker/v/logsdonlab/hgsvc3)](https://hub.docker.com/r/logsdonlab/hgsvc3)

Workflow for [HGSVC3](https://www.internationalgenome.org/human-genome-structural-variation-consortium/) centromere analysis. (WIP)

<img src="docs/rulegraph.svg" width="50%" />

### Usage

#### Clone
```bash
git clone git@github.com:logsdon-lab/hgsvc3.git --recurse-submodules
```

#### Local
```bash
# NOTE: dna-brnn must be installed locally otherwise use --use-singularity
snakemake --use-conda -np --configfile config/config.yaml
```

#### Cluster.
```bash
snakemake -j 100 \
--cluster "bsub -M 40000 -n {threads} -o /dev/null" \
--rerun-triggers mtime \
--configfile config/config.yaml \
--use-conda -n
```

### TODO
* Remove old bedminmax scripts. Only leaving in until output discrepancies finished.
* Add test/move calculate HOR length script to new repo.
