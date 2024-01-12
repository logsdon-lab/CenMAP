# Identify centromeric regions
# Extract centromeric regions


# TODO: Replace with wildcard for assembly
rule intersect_cent_regions:
    input:
        left="GM18989.bed",
        right="cenSat_Annotations_HORs.maxmin.v2.0.500kbp.bed",
    output:
        "GM18989_cens.bed",
    conda:
        "env/env.yaml"
    log:
        "logs/intersect_cent.log",
    shell:
        """
        bedtools intersect -a {input.left} -b {input.right} > {output} 2> {log}
        """


# TODO: This can be parallelized by adding unique fnames with wildcards.
# TODO: What's difference between below?
rule collapse_cens_contigs:
    input:
        regions=rules.intersect_cent_regions.output,
    output:
        "GM18989_1_centromeric_contigs.bed",
    conda:
        "env/env.yaml"
    log:
        "logs/collapse_cent_ctgs.log",
    shell:
        """
        awk -v OFS="\t" '{print $6, $7, $8, $1, $5}' {input.regions} > tmp1.bed

        ./bedminmax.py -i tmp1.bed -o tmp2.bed &> {log}

        awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $3-$2}' tmp2.bed | sort -k1,1 > {output}
        """


# In output dir of asm_to_ref_alignment
# TODO: what's the tail call for?
rule collapse_cens_contigs_only_t2t_cens:
    input:
        regions="GM18989.bed",
    output:
        "GM18989_centromeric_contigs.bed",
    conda:
        "env/env.yaml"
    log:
        "logs/collapse_cent_ctgs_only_t2t_cens.log",
    shell:
        """
        awk -v OFS="\t" '{print $6, $7, $8, $1, $5}' {input.regions} | tail -n +2 > tmp1.bed

        ./bedminmax.py -i tmp1.bed -o tmp2.bed &> {log}

        awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $3-$2}' tmp2.bed | sort -k1,1 > {output}
        """


# TODO: can awk comp be collapsed?
#  awk -F"," -v OFS="," '{dst=($2-$1)} {if (dst > 5 && $2=$3) print $1,$2,$3,$4, dst } ' test.tsv
rule intersect_both_cens_contigs:
    input:
        left="GM18989_centromeric_contigs.bed",
        right="../../../t2t_chm13_v2.0/bed/centromeres/GM18989_centromeric_contigs.bed",
    output:
        "{}_centromeric_regions.all.bed",
    conda:
        "env/env.yaml"
    params:
        # Left-outer-join. Each feat in A, report each overlap w/B.
        # No overlaps, return NULL feat for B
        intersect_params="-loj",
        thr=1000000,
    log:
        "logs/intersect_both_cens_contigs.log",
    shell:
        """
        bedtools intersect {params.intersect_params} \
            -a {input.left} \
            -b {input.right} | \
        awk -v OFS="\t" -F"\t" '{if ($6==$12) print $1, $2, $3, $4, $5, $6, $3-$2}' | \
        awk '$7>{params.thr}' | \
        uniq > tmp3.bed 2> {log}

        # Collapse ctgs grouped by label.
        ./bedminmax.py -i tmp3.bed -o tmp4.bed 2> {log}

        bedtools intersect {params.intersect_params} -a tmp4.bed -b tmp3.bed | \
        awk -v OFS="\t" '{print $1, $2, $3, $4, $9, $10, $3-$2}' | \
        sort -k5,5 | \
        uniq > {output} 2> {log}
        """


# Reorient the regions.
rule extract_fwd_rev_regions:
    input:
        all_regions="GM18989_centromeric_regions.all.bed",
        combined_assembly="/net/eichler/vol27/projects/hgsvc/nobackups/analysis/glennis/assemblies/GM18989.vrk-ps-sseq.asm-combined.fa",
    output:
        fwd_cen_regions="GM18989_centromeric_regions.fwd.bed",
        rev_cen_regions="GM18989_centromeric_regions.rev.bed",
        fwd_cen_seq="GM18989_centromeric_regions.fwd.fa",
        rev_cen_seq="GM18989_centromeric_regions.rev.fa",
    log:
        "logs/extract_fwd_rev_regions.log",
    conda:
        "env/env.yaml"
    shell:
        """
        awk -v OFS="\t" '{if ($6=="+") print}' {input.all_regions} > {output.fwd_cen_regions}
        awk -v OFS="\t" '{if ($6=="-") print}' {input.all_regions} > {output.rev_cen_regions}
        seqtk subseq {input.combined_assembly} {output.fwd_cen_regions} > {output.fwd_cen_seq}
        seqtk subseq {input.combined_assembly} {output.rev_cen_regions} | seqtk seq -r - > {output.rev_cen_seq}
        """


# TODO: Move to start before alignment.
rule create_fwd_ctg_name_legend:
    input:
        regions=rules.extract_fwd_rev_regions.fwd_cen_regions,
    output:
        "{}.legend.fwd.txt",
    log:
        "logs/create_fwd_ctg_name_legend.log",
    conda:
        "env/env.yaml"
    shell:
        """
        awk -v OFS="\t" '{print $0, FILENAME}' {input.regions} | \
        sed 's/_/\t/g' | \
        awk -v OFS="\t" '{print $5, $8"_"$4"_"$5}' | \
        sort -k2,2 > {output} 2> {log}
        """


use rule create_fwd_ctg_name_legend as create_rev_ctg_name_legend with:
    input:
        regions=rules.extract_fwd_rev_regions.rev_cen_regions,
    output:
        "{}.legend.rev.txt",
    log:
        "logs/create_rev_ctg_name_legend.log",


rule split_fwd_cens_assembly_fasta:
    input:
        rules.extract_fwd_rev_regions.fwd_cen_seq,
    output:
        "{}.fwd.txt",
    conda:
        "env/env.yaml"
    log:
        "logs/split_fwd_cens_assembly_fasta.log",
    shell:
        """
        sed 's/>/>\n/g' {input} | sed 's/:/\n/g' > {output}
        """


use rule split_fwd_cens_assembly_fasta as split_rev_cens_assembly_fasta with:
    input:
        rules.extract_fwd_rev_regions.rev_cen_seq,
    output:
        "{}.rev.txt",
    log:
        "logs/split_rev_cens_assembly_fasta.log",


rule rename_cens_fwd_ctgs:
    input:
        # a
        legend=rules.create_fwd_ctg_name_legend,
        # b
        seq=rules.split_fwd_cens_assembly_fasta,
    output:
        "{}_centromeric_regions.renamed.fwd.fa",
    conda:
        "env/env.yaml"
    log:
        "logs/rename_cens_fwd_ctgs.log",
    shell:
        """
        awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2;next} {print ($1 in a ? a[$1] : $1)}' {input.legend} {input.seq} | \
        awk '{printf "%s%s", (/>/ ? ors : OFS), $0; ors=ORS} END{print ":"}' | \
        sed 's/> />/g' | \
        sed 's/\([0-9]\) \([0-9]\)/\1:\2/g' | \
        tr " " "\n" > {output} 2> {log}
        """


use rule rename_cens_fwd_ctgs as rename_cens_rev_ctgs with:
    input:
        legend=rules.create_rev_ctg_name_legend,
        seq=rules.split_rev_cens_assembly_fasta,
    output:
        "{}_centromeric_regions.renamed.rev.fa",
    log:
        "logs/rename_cens_rev_ctgs.log",


rule index_renamed_cens_ctgs:
    input:
        rules.rename_cens_fwd_ctgs.output,
    output:
        "{}_centromeric_regions.renamed.rev.fai",
    conda:
        "env/env.yaml"
    shell:
        """
        samtools faidx {input}
        """
