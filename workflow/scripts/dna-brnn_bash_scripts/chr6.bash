#!/bin/env bash

rm chr6_tmp.*.bed

grep "chr6:" chm13_cens.trimmed.bed | sed 's/:/\t/g' | sed 's/-/\t/g' | awk -v OFS="\t" '{print $1, $2+$4, $2+$5, $6, $5-$4}' | awk '$4==2' | awk '$5>1000' > chr6_tmp.fwd.bed

for sample in $(ls ./ | grep "bed" | grep "centromeric_regions" | awk -F "." '{print $1}' | sort | uniq); do
	grep "chr6_" ${sample}.renamed.fwd.bed | sed 's/:/\t/g' | sed 's/-/\t/g' | awk -v OFS="\t" '{print $1"-"$2, $3+$5, $3+$6, $7, $6-$5}' | awk '$4==2' | awk '$5>1000' >> chr6_tmp.fwd.bed
	/net/eichler/vol28/home/glogsdon/utilities/bedminmax.py -i chr6_tmp.fwd.bed | awk -v OFS="\t" '{print $1, $2, $3, $3-$2}' | awk -v OFS="\t" '{print $1, $2-493250, $3+434399, $3-$2}' | awk '$4>1000000' | awk -v OFS="\t" '$2<0 {$2=0}1' > chr6_contigs.fwd.ALR.bed
	grep "chr6_" ${sample}.renamed.rev.bed | sed 's/:/\t/g' | sed 's/-/\t/g' | awk -v OFS="\t" '{print $1"-"$2, $4-$6, $4-$5, $7, $6-$5}' | awk '$4==2' | awk '$5>1000' >> chr6_tmp.rev.bed
	/net/eichler/vol28/home/glogsdon/utilities/bedminmax.py -i chr6_tmp.rev.bed | awk -v OFS="\t" '{print $1, $2, $3, $3-$2}' | awk -v OFS="\t" '{print $1, $2-493250, $3+434399, $3-$2}' | awk '$4>1000000' | awk -v OFS="\t" '$2<0 {$2=0}1' > chr6_contigs.rev.ALR.bed
done

rm chr6_tmp.*.bed