#!/bin/env bash

rm chr11_tmp.*.bed

grep "chr11:" chm13_cens.trimmed.bed | sed 's/:/\t/g' | sed 's/-/\t/g' | awk -v OFS="\t" '{print $1, $2+$4, $2+$5, $6, $5-$4}' | awk '$4==2' | awk '$5>1000' > chr11_tmp.fwd.bed

for sample in $(ls ./ | grep "bed" | grep "centromeric_regions" | awk -F "." '{print $1}' | sort | uniq); do
	grep "chr11_" ${sample}.renamed.fwd.bed | sed 's/:/\t/g' | sed 's/-/\t/g' | awk -v OFS="\t" '{print $1"-"$2, $3+$5, $3+$6, $7, $6-$5}' | awk '$4==2' | awk '$5>1000' >> chr11_tmp.fwd.bed
	/net/eichler/vol28/home/glogsdon/utilities/bedminmax.py -i chr11_tmp.fwd.bed | awk -v OFS="\t" '{print $1, $2, $3, $3-$2}' | awk -v OFS="\t" '{print $1, $2-519650, $3+502250, $3-$2}' | awk '$4>1000000' | awk -v OFS="\t" '$2<0 {$2=0}1' > chr11_contigs.fwd.ALR.bed
	grep "chr11_" ${sample}.renamed.rev.bed | sed 's/:/\t/g' | sed 's/-/\t/g' | awk -v OFS="\t" '{print $1"-"$2, $4-$6, $4-$5, $7, $6-$5}' | awk '$4==2' | awk '$5>1000' >> chr11_tmp.rev.bed
	/net/eichler/vol28/home/glogsdon/utilities/bedminmax.py -i chr11_tmp.rev.bed | awk -v OFS="\t" '{print $1, $2, $3, $3-$2}' | awk -v OFS="\t" '{print $1, $2-519650, $3+502250, $3-$2}' | awk '$4>1000000' | awk -v OFS="\t" '$2<0 {$2=0}1' > chr11_contigs.rev.ALR.bed
done

rm chr11_tmp.*.bed
