#!/bin/bash

set -euo pipefail

ctg_name="HG00731_chrY_h2tg000043l#1-12137861"

samtools view results_hifiasm/nucflag/HG00731_hifi.bam $ctg_name | cut -f 1 | sort | uniq > ${ctg_name}_reads.list

for read_file in data/raw_data_link/HG00731/*.fastq.gz; do
    seqtk subseq $read_file ${ctg_name}_reads.list >> $ctg_name.fastq
done

seqkit rmdup $ctg_name.fastq | bgzip > ${ctg_name}_reads.fastq.gz
