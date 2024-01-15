#!/bin/bash

set -euo pipefail

ref=data/chm13_t2t/nobackups/assemblies/chm13_v1.1_plus38Y_par_masked.fasta
asat_annot=data/AlphaSatelliteMapping/nobackups/FindingAlphaSat/t2t_chr8/ape_cens/chm1/asat_annotation/chm13_v2.0_hor_arrays.bed
cens_regs=data/AlphaSatelliteMapping/nobackups/FindingAlphaSat/t2t_chr8/chm13/centromeres/annotations/cenSat_Annotations_HORs.maxmin.v2.0.500kbp.bed
mkdir -p $(dirname $ref) && touch $ref
mkdir -p $(dirname $asat_annot) && touch $asat_annot
mkdir -p $(dirname $cens_regs) && touch $cens_regs

# Make fake assemblies
touch data/HG002_1.fa.gz
touch data/HG002_2.fa.gz
touch data/HG00733_1.fa.gz
touch data/HG00733_2.fa.gz

# Make fake combined verkko assemblies
assembly_dir=data/assemblies
mkdir -p $assembly_dir
awk 'NR > 1 { print $1 }' config/table.asm.tbl | \
    xargs -I [] \
    bash -c "
        touch $assembly_dir/[]_1.vrk-ps-sseq.asm-combined.fa && \
        touch $assembly_dir/[]_2.vrk-ps-sseq.asm-combined.fa
    "

touch data/AS-HORs-hmmer3.0-170921.hmm
