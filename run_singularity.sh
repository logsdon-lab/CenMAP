#!/bin/bash

set -euo pipefail

singularity run hgsvc3_latest.sif \
--use-conda -p \
--rerun-triggers mtime \
--configfile config/config.yaml \
--notemp "$@"
