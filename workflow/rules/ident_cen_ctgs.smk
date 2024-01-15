# Identify centromeric regions
# Extract centromeric regions


# Get centromeric regions from alignments to t2t-chm13
rule intersect_cen_regions:
    input:
        left="{sample}.bed",
        right=config["align_asm_to_ref"]["cens_500kbp_regions"],
    output:
        cen_regions="{sample}_cens.bed",
    conda:
        "../env/tools.yaml"
    log:
        "logs/intersect_cent_{sample}.log",
    shell:
        """
        bedtools intersect -a {input.left} -b {input.right} > {output} 2> {log}
        """


# TODO: This can be parallelized by adding unique fnames with wildcards.
# TODO: What's difference between below?
rule collapse_cens_contigs:
    input:
        script="workflow/scripts/bedmin_max.py",
        regions=rules.intersect_cen_regions.output,
    output:
        "{sample}_1_centromeric_contigs.bed",
    conda:
        "../env/py.yaml"
    log:
        "logs/collapse_cent_ctgs_{sample}.log",
    shell:
        """
        awk -v OFS="\t" '{print $6, $7, $8, $1, $5}' {input.regions} > tmp1.bed

        python {input.script} -i tmp1.bed -o tmp2.bed &> {log}

        awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $3-$2}' tmp2.bed | sort -k1,1 > {output}
        """


# In output dir of asm_to_ref_alignment
# TODO: what's the tail call for?
rule collapse_cens_contigs_only_t2t_cens:
    input:
        script="workflow/scripts/bedmin_max.py",
        regions="{sample}.bed",
    output:
        "{sample}_centromeric_contigs.bed",
    conda:
        "../env/py.yaml"
    log:
        "logs/collapse_cent_ctgs_only_t2t_cens_{sample}.log",
    shell:
        """
        awk -v OFS="\t" '{print $6, $7, $8, $1, $5}' {input.regions} | tail -n +2 > tmp1.bed

        python {input.script} -i tmp1.bed -o tmp2.bed &> {log}

        awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $3-$2}' tmp2.bed | sort -k1,1 > {output}
        """


# TODO: can awk comp be collapsed?
#  awk -F"," -v OFS="," '{dst=($2-$1)} {if (dst > 5 && $2=$3) print $1,$2,$3,$4, dst } ' test.tsv
rule intersect_both_cens_contigs:
    input:
        left="{sample}_centromeric_contigs.bed",
        right="../../../t2t_chm13_v2.0/bed/centromeres/GM18989_centromeric_contigs.bed",
    output:
        "{sample}_centromeric_regions.all.bed",
    conda:
        "../env/env.yaml"
    params:
        # Left-outer-join. Each feat in A, report each overlap w/B.
        # No overlaps, return NULL feat for B
        intersect_params="-loj",
        thr=1_000_00,
    log:
        "logs/intersect_both_cens_contigs_{sample}.log",
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
        all_regions="{sample}_centromeric_regions.all.bed",
        combined_assembly="/net/eichler/vol27/projects/hgsvc/nobackups/analysis/glennis/assemblies/GM18989.vrk-ps-sseq.asm-combined.fa",
    output:
        fwd_cen_regions="{sample}_centromeric_regions.fwd.bed",
        rev_cen_regions="{sample}_centromeric_regions.rev.bed",
        fwd_cen_seq="{sample}_centromeric_regions.fwd.fa",
        rev_cen_seq="{sample}_centromeric_regions.rev.fa",
    log:
        "logs/extract_fwd_rev_regions.log",
    conda:
        "../env/tools.yaml"
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
        regions=rules.extract_fwd_rev_regions.output.fwd_cen_regions,
    output:
        "{}.legend.fwd.txt",
    log:
        "logs/create_fwd_ctg_name_legend.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        awk -v OFS="\t" '{print $0, FILENAME}' {input.regions} | \
        sed 's/_/\t/g' | \
        awk -v OFS="\t" '{print $5, $8"_"$4"_"$5}' | \
        sort -k2,2 > {output} 2> {log}
        """


use rule create_fwd_ctg_name_legend as create_rev_ctg_name_legend with:
    input:
        regions=rules.extract_fwd_rev_regions.output.rev_cen_regions,
    output:
        "{}.legend.rev.txt",
    log:
        "logs/create_rev_ctg_name_legend.log",


rule split_fwd_cens_assembly_fasta:
    input:
        rules.extract_fwd_rev_regions.output.fwd_cen_seq,
    output:
        "{}.fwd.txt",
    conda:
        "../env/tools.yaml"
    log:
        "logs/split_fwd_cens_assembly_fasta.log",
    shell:
        """
        sed 's/>/>\n/g' {input} | sed 's/:/\n/g' > {output}
        """


use rule split_fwd_cens_assembly_fasta as split_rev_cens_assembly_fasta with:
    input:
        rules.extract_fwd_rev_regions.output.rev_cen_seq,
    output:
        "{}.rev.txt",
    log:
        "logs/split_rev_cens_assembly_fasta.log",


rule rename_cens_fwd_ctgs:
    input:
        # a
        legend=rules.create_fwd_ctg_name_legend.output,
        # b
        seq=rules.split_fwd_cens_assembly_fasta.output,
    output:
        "{}_centromeric_regions.renamed.fwd.fa",
    conda:
        "../env/tools.yaml"
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
        legend=rules.create_rev_ctg_name_legend.output,
        seq=rules.split_rev_cens_assembly_fasta.output,
    output:
        "{}_centromeric_regions.renamed.rev.fa",
    log:
        "logs/rename_cens_rev_ctgs.log",


rule index_renamed_cens_ctgs:
    input:
        rules.rename_cens_fwd_ctgs.output,
    output:
        "{sample}_centromeric_regions.renamed.rev.fai",
    conda:
        "../env/tools.yaml"
    shell:
        """
        samtools faidx {input}
        """


rule ident_cen_ctgs_all:
    input:
        expand(),
