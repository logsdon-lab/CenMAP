#!/bin/bash

set -euo pipefail

ref=data/reference/T2T-CHM13v2.fasta
mkdir -p $(dirname $ref) && touch $ref

# Make fake assemblies
touch data/HG00171_1.fa.gz
touch data/HG00171_2.fa.gz

# Make fake combined verkko assemblies
assembly_dir=data/assemblies
mkdir -p $assembly_dir
awk 'NR > 1 { print $1 }' config/table.asm.tbl | \
    xargs -I [] \
    bash -c "
        mkdir -p $assembly_dir/[]
        touch $assembly_dir/[]/[].vrk-ps-sseq.asm-hap1.fasta.gz && \
        touch $assembly_dir/[]/[].vrk-ps-sseq.asm-hap2.fasta.gz
    "

touch data/AS-HORs-hmmer3.0-170921.hmm
