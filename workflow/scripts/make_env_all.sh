#!/bin/bash

set -euo pipefail

conda-merge $(snakemake -p -c 1 --sdm apptainer --configfile test/config/config.yaml --containerize | grep COPY | awk '{ print $2}') | \
    grep -v "name:" > env_all_.yaml
