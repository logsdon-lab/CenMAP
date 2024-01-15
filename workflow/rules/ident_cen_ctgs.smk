# Identify centromeric regions
# Extract centromeric regions

REF_NAME = get_ref_name()
OUTPUT_DIR = config["ident_cen_ctgs"]["output_dir"]


# Get centromeric regions from alignments to t2t-chm13 ONLY
rule intersect_cen_regions:
    input:
        left=lambda wc: expand(
            rules.ref_align_aln_to_bed.output, ref=[REF_NAME], sm=wc.sm
        ),
        right=config["align_asm_to_ref"]["cens_500kbp_regions"],
    output:
        cen_regions=os.path.join(f"results/{REF_NAME}/bed", "{sm}_cens.bed"),
    conda:
        "../env/tools.yaml"
    log:
        "logs/intersect_cent_{sm}.log",
    shell:
        """
        bedtools intersect -a {input.left} -b {input.right} > {output} 2> {log}
        """


rule collapse_cens_contigs:
    input:
        script="workflow/scripts/bedmin_max.py",
        regions=rules.intersect_cen_regions.output,
    output:
        os.path.join(f"results/{REF_NAME}/bed", "{sm}_centromeric_contigs.bed"),
    conda:
        "../env/py.yaml"
    log:
        "logs/collapse_cent_ctgs_{sm}.log",
    shell:
        """
        awk -v OFS="\t" '{{print $6, $7, $8, $1, $5}}' {input.regions} > tmp1.bed

        python {input.script} -i tmp1.bed -o tmp2.bed &> {log}

        awk -v OFS="\t" '{{print $1, $2, $3, $4, $5, $3-$2}}' tmp2.bed | sort -k1,1 > {output}
        """


# T2T-CHM13 CENS ONLY
# In output dir of asm_to_ref_alignment
# TODO: what's the tail call for?
rule collapse_cens_contigs_only_t2t_cens:
    input:
        script="workflow/scripts/bedmin_max.py",
        regions=lambda wc: expand(
            rules.ref_align_aln_to_bed.output, ref=[REF_NAME], sm=wc.sm
        ),
    output:
        os.path.join(f"results/{REF_NAME}_cens/bed", "{sm}_centromeric_contigs.bed"),
    conda:
        "../env/py.yaml"
    log:
        "logs/collapse_cent_ctgs_only_t2t_cens_{sm}.log",
    shell:
        """
        awk -v OFS="\t" '{{print $6, $7, $8, $1, $5}}' {input.regions} | tail -n +2 > tmp1.bed

        python {input.script} -i tmp1.bed -o tmp2.bed &> {log}

        awk -v OFS="\t" '{{print $1, $2, $3, $4, $5, $3-$2}}' tmp2.bed | sort -k1,1 > {output}
        """


# TODO: can awk comp be collapsed?
#  awk -F"," -v OFS="," '{dst=($2-$1)} {if (dst > 5 && $2=$3) print $1,$2,$3,$4, dst } ' test.tsv
rule intersect_filter_both_cens_ctgs:
    input:
        left=f"results/{REF_NAME}_cens/bed/{sm}_centromeric_contigs.bed",
        # TODO: Ask about 1
        right=f"results/{REF_NAME}/bed/{sm}_centromeric_contigs.bed",
    output:
        filt_regions=temp(
            os.path.join(OUTPUT_DIR, "{sm}_centromeric_regions.inter_filt.bed")
        ),
    conda:
        "../env/tools.yaml"
    params:
        # Left-outer-join. Each feat in A, report each overlap w/B.
        # No overlaps, return NULL feat for B
        intersect_params="-loj",
        thr=1_000_00,
    log:
        "logs/intersect_filter_both_cens_ctgs_{sm}.log",
    shell:
        """
        bedtools intersect {params.intersect_params} \
            -a {input.left} \
            -b {input.right} | \
        awk -v OFS="\t" -F"\t" '{{if ($6==$12) print $1, $2, $3, $4, $5, $6, $3-$2}}' | \
        awk '$7>{params.thr}' | \
        uniq > {output} 2> {log}
        """


rule collapse_intersected_filtered_cen_ctgs:
    input:
        script="workflow/scripts/bedmin_max.py",
        filt_regions=rules.intersect_filter_both_cens_ctgs.output,
    output:
        collapsed_regions=temp(
            os.path.join(OUTPUT_DIR, "{sm}_centromeric_regions.collapsed.bed")
        ),
    log:
        "logs/collapse_intersected_filtered_cen_ctgs_{sm}.log",
    conda:
        "../env/py.yaml"
    shell:
        """
        # Collapse ctgs grouped by label.
        python {input.script} -i {input.filt_regions} -o {output.collapsed_regions} 2> {log}
        """


rule reintersect_sort_uniq_cens_ctgs:
    input:
        collapsed_regions=rules.collapse_intersected_filtered_cen_ctgs.output,
        filt_regions=rules.intersect_filter_both_cens_ctgs.output,
    output:
        os.path.join(OUTPUT_DIR, "{sm}_centromeric_regions.all.bed"),
    params:
        intersect_params="-loj",
    conda:
        "../env/tools.yaml"
    log:
        "logs/reintersect_sort_uniq_cens_ctgs_{sm}.log",
    shell:
        """
        bedtools intersect {params.intersect_params} -a {input.collapsed_regions} -b {input.filt_regions}| \
        awk -v OFS="\t" '{{print $1, $2, $3, $4, $9, $10, $3-$2}}' | \
        sort -k5,5 | \
        uniq > {output} 2> {log}
        """


# Reorient the regions.
rule extract_fwd_rev_regions:
    input:
        all_regions=rules.reintersect_sort_uniq_cens_ctgs.output,
        combined_assembly=os.path.join(
            config["ident_cen_ctgs"]["comb_assemblies_dir"],
            "{sm}.vrk-ps-sseq.asm-combined.fa",
        ),
        # combined_assembly="/net/eichler/vol27/projects/hgsvc/nobackups/analysis/glennis/assemblies/{sm}.vrk-ps-sseq.asm-combined.fa",
    output:
        fwd_cen_regions=os.path.join(OUTPUT_DIR, "{sm}_centromeric_regions.fwd.bed"),
        rev_cen_regions=os.path.join(OUTPUT_DIR, "{sm}_centromeric_regions.rev.bed"),
        fwd_cen_seq=os.path.join(OUTPUT_DIR, "{sm}_centromeric_regions.fwd.fa"),
        rev_cen_seq=os.path.join(OUTPUT_DIR, "{sm}_centromeric_regions.rev.fa"),
    log:
        "logs/extract_fwd_rev_regions_{sm}.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        awk -v OFS="\t" '{{if ($6=="+") print}}' {input.all_regions} > {output.fwd_cen_regions}
        awk -v OFS="\t" '{{if ($6=="-") print}}' {input.all_regions} > {output.rev_cen_regions}
        seqtk subseq {input.combined_assembly} {output.fwd_cen_regions} > {output.fwd_cen_seq}
        seqtk subseq {input.combined_assembly} {output.rev_cen_regions} | seqtk seq -r - > {output.rev_cen_seq}
        """


# TODO: Move to start before alignment.
rule create_fwd_ctg_name_legend:
    input:
        regions=rules.extract_fwd_rev_regions.output.fwd_cen_regions,
    output:
        os.path.join(OUTPUT_DIR, "{sm}.legend.fwd.txt"),
    log:
        "logs/create_fwd_ctg_name_legend_{sm}.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        awk -v OFS="\t" '{{print $0, FILENAME}}' {input.regions} | \
        sed 's/_/\t/g' | \
        awk -v OFS="\t" '{{print $5, $8"_"$4"_"$5}}' | \
        sort -k2,2 > {output} 2> {log}
        """


use rule create_fwd_ctg_name_legend as create_rev_ctg_name_legend with:
    input:
        regions=rules.extract_fwd_rev_regions.output.rev_cen_regions,
    output:
        os.path.join(OUTPUT_DIR, "{sm}.legend.rev.txt"),
    log:
        "logs/create_rev_ctg_name_legend_{sm}.log",


rule split_fwd_cens_assembly_fasta:
    input:
        rules.extract_fwd_rev_regions.output.fwd_cen_seq,
    output:
        os.path.join(OUTPUT_DIR, "{sm}.fwd.txt"),
    conda:
        "../env/tools.yaml"
    log:
        "logs/split_fwd_cens_assembly_fasta_{sm}.log",
    shell:
        """
        sed 's/>/>\n/g' {input} | sed 's/:/\n/g' > {output}
        """


use rule split_fwd_cens_assembly_fasta as split_rev_cens_assembly_fasta with:
    input:
        rules.extract_fwd_rev_regions.output.rev_cen_seq,
    output:
        os.path.join(OUTPUT_DIR, "{sm}.rev.txt"),
    log:
        "logs/split_rev_cens_assembly_fasta_{sm}.log",


rule rename_cens_fwd_ctgs:
    input:
        # a
        legend=rules.create_fwd_ctg_name_legend.output,
        # b
        seq=rules.split_fwd_cens_assembly_fasta.output,
    output:
        os.path.join(OUTPUT_DIR, "{sm}_centromeric_regions.renamed.fwd.fa"),
    conda:
        "../env/tools.yaml"
    log:
        "logs/rename_cens_fwd_ctgs_{sm}.log",
    shell:
        """
        awk 'BEGIN{{FS=OFS="\t"}} NR==FNR{{a[$1]=$2;next}} {{print ($1 in a ? a[$1] : $1)}}' {input.legend} {input.seq} | \
        awk '{{printf "%s%s", (/>/ ? ors : OFS), $0; ors=ORS}} END{{print ":"}}' | \
        sed 's/> />/g' | \
        sed 's/\([0-9]\) \([0-9]\)/\1:\2/g' | \
        tr " " "\n" > {output} 2> {log}
        """


use rule rename_cens_fwd_ctgs as rename_cens_rev_ctgs with:
    input:
        legend=rules.create_rev_ctg_name_legend.output,
        seq=rules.split_rev_cens_assembly_fasta.output,
    output:
        os.path.join(OUTPUT_DIR, "{sm}_centromeric_regions.renamed.rev.fa"),
    log:
        "logs/rename_cens_rev_ctgs_{sm}.log",


rule index_renamed_cens_ctgs:
    input:
        fwd=rules.rename_cens_fwd_ctgs.output,
        rev=rules.rename_cens_rev_ctgs.output,
    output:
        fwd=os.path.join(OUTPUT_DIR, "{sm}_centromeric_regions.renamed.fwd.fai"),
        rev=os.path.join(OUTPUT_DIR, "{sm}_centromeric_regions.renamed.rev.fai"),
    conda:
        "../env/tools.yaml"
    log:
        "logs/index_renamed_cens_ctgs_{sm}.log",
    shell:
        """
        samtools faidx {input} 2> {log}
        """


rule ident_cen_ctgs_all:
    input:
        expand(rules.intersect_cen_regions.output, sm=SAMPLES.index),
        expand(rules.collapse_cens_contigs.output, sm=SAMPLES.index),
        expand(
            rules.collapse_cens_contigs_only_t2t_cens.output,
            sm=SAMPLES.index,
        ),
        expand(rules.intersect_filter_both_cens_ctgs.output, sm=SAMPLES.index),
        expand(rules.collapse_intersected_filtered_cen_ctgs.output, sm=SAMPLES.index),
        expand(rules.reintersect_sort_uniq_cens_ctgs.output, sm=SAMPLES.index),
        expand(rules.extract_fwd_rev_regions.output, sm=SAMPLES.index),
        # Renaming section
        expand(rules.create_fwd_ctg_name_legend.output, sm=SAMPLES.index),
        expand(rules.create_rev_ctg_name_legend.output, sm=SAMPLES.index),
        expand(rules.split_fwd_cens_assembly_fasta.output, sm=SAMPLES.index),
        expand(rules.split_rev_cens_assembly_fasta.output, sm=SAMPLES.index),
        expand(rules.rename_cens_fwd_ctgs.output, sm=SAMPLES.index),
        expand(rules.rename_cens_rev_ctgs.output, sm=SAMPLES.index),
        expand(rules.index_renamed_cens_ctgs.output, sm=SAMPLES.index),
