#!/bin/env bash

rm chr8_tmp.*.bed

grep "chr8:" chm13_cens.trimmed.bed | sed 's/:/\t/g' | sed 's/-/\t/g' | awk -v OFS="\t" '{print $1, $2+$4, $2+$5, $6, $5-$4}' | awk '$4==2' | awk '$5>1000' > chr8_tmp.fwd.bed

for sample in $(ls ./ | grep "bed" | grep "centromeric_regions" | awk -F "." '{print $1}' | sort | uniq); do
	grep "chr8_" ${sample}.renamed.fwd.bed | sed 's/:/\t/g' | sed 's/-/\t/g' | awk -v OFS="\t" '{print $1"-"$2, $3+$5, $3+$6, $7, $6-$5}' | awk '$4==2' | awk '$5>1000' >> chr8_tmp.fwd.bed
	/net/eichler/vol28/home/glogsdon/utilities/bedminmax.py -i chr8_tmp.fwd.bed | awk -v OFS="\t" '{print $1, $2, $3, $3-$2}' | awk -v OFS="\t" '{print $1, $2-496550, $3+529550, $3-$2}' | awk '$4>1000000' | awk -v OFS="\t" '$2<0 {$2=0}1' > chr8_contigs.fwd.ALR.bed
	grep "chr8_" ${sample}.renamed.rev.bed | sed 's/:/\t/g' | sed 's/-/\t/g' | awk -v OFS="\t" '{print $1"-"$2, $4-$6, $4-$5, $7, $6-$5}' | awk '$4==2' | awk '$5>1000' >> chr8_tmp.rev.bed
	/net/eichler/vol28/home/glogsdon/utilities/bedminmax.py -i chr8_tmp.rev.bed | awk -v OFS="\t" '{print $1, $2, $3, $3-$2}' | awk -v OFS="\t" '{print $1, $2-496550, $3+529550, $3-$2}' | awk '$4>1000000' | awk -v OFS="\t" '$2<0 {$2=0}1' > chr8_contigs.rev.ALR.bed
done

rm chr8_tmp.*.bed
