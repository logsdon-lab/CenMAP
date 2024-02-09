#!/bin/bash

set -euo pipefail

num_cores=12
max_mem=61440

bsub \
-e test_hgsvc3.err.log -o test_hgsvc3.out.log \
-n $num_cores \
-M $max_mem \
-R "rusage [mem=$max_mem] span[hosts=1]" \
singularity run hgsvc3_latest.sif \
--use-conda -p \
-c $num_cores \
--rerun-triggers mtime \
--configfile config/config.yaml \
--notemp "$@"
