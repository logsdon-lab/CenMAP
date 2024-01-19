# Identify centromeric regions
# Extract centromeric regions

OUTPUT_DIR_T2T_REF_REGIONS = f"results/{REF_NAME}/bed"
OUTPUT_DIR_T2T_REF_CENS_ONLY_REGIONS = f"results/{REF_NAME}_cens/bed"


# Get centromeric regions from alignments to t2t-chm13 ONLY
rule intersect_cen_regions:
    input:
        left=lambda wc: expand(
            rules.ref_align_aln_to_bed.output, ref=[REF_NAME], sm=wc.sm
        ),
        right=config["align_asm_to_ref"]["cens_500kbp_regions"],
    output:
        cen_regions=os.path.join(OUTPUT_DIR_T2T_REF_REGIONS, "{sm}_cens.bed"),
    conda:
        "../env/tools.yaml"
    log:
        "logs/intersect_cent_{sm}.log",
    shell:
        """
        bedtools intersect -a {input.left} -b {input.right} > {output} 2> {log}
        """


# 6. query_name
# 7. query_start
# 8. query_end
# 4. reference_length
# 1. reference_name
# 5. strand
rule collapse_cens_contigs:
    input:
        # TODO: This one has orientation, v2 doesn't. What's up with that?
        script="workflow/scripts/bedminmax.py",
        regions=rules.intersect_cen_regions.output,
    output:
        refmt_regions=temp(
            os.path.join(
                OUTPUT_DIR_T2T_REF_REGIONS, "len_calc_{sm}_centromeric_contigs.bed"
            )
        ),
        regions=os.path.join(OUTPUT_DIR_T2T_REF_REGIONS, "{sm}_centromeric_contigs.bed"),
    conda:
        "../env/py.yaml"
    params:
        # NR>1 if yes, blank if no.
        awk_print_cmd="'{ print $6, $7, $8, $4, $1, $5 }'",
    log:
        "logs/collapse_cent_ctgs_{sm}.log",
    shell:
        """
        awk -v OFS="\\t" {params.awk_print_cmd} \
        {input.regions} > {output.refmt_regions} 2> {log}

        # Calculate length of ref region.
        {{ python {input.script} -i {output.refmt_regions} | \
        awk -v OFS="\\t" '{{print $0, $3-$2}}' | \
        sort -k1,1;}} > {output.regions} 2> {log}
        """


# 1st awk
# 6. query_name
# 7. query_start
# 8. query_end
# 9. query_length
# 1. reference_name
# 5. strand


# 2nd awk (# relative to 1st)
# 1. query_name
# 2. query_start
# 3. query_end
# 4. query_length
# 5. reference_name
# 6. strand
# +. sub(query_end, query_start)
# T2T-CHM13 CENS ONLY
use rule collapse_cens_contigs as collapse_cens_contigs_only_t2t_cens with:
    input:
        script="workflow/scripts/bedminmax.py",
        # asm aligned to t2t_chm13_v2.hor_arrays_masked.500kbp.fa
        regions=lambda wc: expand(
            rules.ref_align_aln_to_bed.output, ref=[f"{REF_NAME}_cens"], sm=wc.sm
        ),
    output:
        refmt_regions=temp(
            os.path.join(
                OUTPUT_DIR_T2T_REF_CENS_ONLY_REGIONS,
                "len_calc_{sm}_centromeric_contigs.bed",
            )
        ),
        regions=os.path.join(
            OUTPUT_DIR_T2T_REF_CENS_ONLY_REGIONS, "{sm}_centromeric_contigs.bed"
        ),
    params:
        awk_print_cmd="'NR>1 {print $6, $7, $8, $9, $1, $5}'",
    log:
        "logs/collapse_cent_ctgs_only_t2t_cens_{sm}.log",


# TODO: can awk comp be collapsed?
#  awk -F"," -v OFS="," '{dst=($2-$1)} {if (dst > 5 && $2=$3) print $1,$2,$3,$4, dst } ' test.tsv
rule intersect_filter_both_cens_ctgs:
    input:
        left=os.path.join(
            OUTPUT_DIR_T2T_REF_CENS_ONLY_REGIONS, "{sm}_centromeric_contigs.bed"
        ),
        right=os.path.join(OUTPUT_DIR_T2T_REF_REGIONS, "{sm}_centromeric_contigs.bed"),
    output:
        filt_regions=temp(
            os.path.join(
                config["ident_cen_ctgs"]["output_dir"],
                "{sm}_centromeric_regions.inter_filt.bed",
            )
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
        {{ bedtools intersect {params.intersect_params} \
            -a {input.left} \
            -b {input.right} | \
        awk -v OFS="\\t" -F"\\t" '{{if ($6==$13) print $1, $2, $3, $4, $5, $6, $3-$2}}' | \
        awk '$7>{params.thr}' | \
        uniq;}} > {output} 2> {log}
        """


rule collapse_intersected_filtered_cen_ctgs:
    input:
        script="workflow/scripts/bed-er-minmax.py",
        filt_regions=rules.intersect_filter_both_cens_ctgs.output,
    output:
        collapsed_regions=temp(
            os.path.join(
                config["ident_cen_ctgs"]["output_dir"],
                "{sm}_centromeric_regions.collapsed.bed",
            )
        ),
    params:
        input_cols=" ".join(
            [
                "chr",
                "start",
                "end",
                "length",
                "name",
                "orientation",
                "query_start_end_diff",
            ]
        ),
    log:
        "logs/collapse_intersected_filtered_cen_ctgs_{sm}.log",
    conda:
        "../env/py.yaml"
    shell:
        """
        # Collapse ctgs grouped by label.
        python {input.script} \
        -i {input.filt_regions} \
        -o {output.collapsed_regions} \
        -ci {params.input_cols} 2> {log}
        """


# Exp:
# coll_name
# coll_start
# coll_end
# coll_length
# inter_name
# inter_strand
# coll_diff
rule reintersect_sort_uniq_cens_ctgs:
    input:
        collapsed_regions=rules.collapse_intersected_filtered_cen_ctgs.output,
        filt_regions=rules.intersect_filter_both_cens_ctgs.output,
    output:
        os.path.join(
            config["ident_cen_ctgs"]["output_dir"], "{sm}_centromeric_regions.all.bed"
        ),
    params:
        intersect_params="-loj",
    conda:
        "../env/tools.yaml"
    log:
        "logs/reintersect_sort_uniq_cens_ctgs_{sm}.log",
    shell:
        """
        {{ bedtools intersect {params.intersect_params} -a {input.collapsed_regions} -b {input.filt_regions}| \
        awk -v OFS="\\t" '{{print $1, $2, $3, $4, $11, $12, $3-$2}}' | \
        sort -k5,5 | \
        uniq;}} > {output} 2> {log}
        """


# Reorient the regions.
rule extract_fwd_rev_regions:
    input:
        all_regions=rules.reintersect_sort_uniq_cens_ctgs.output,
        # Weird hack because sm wildcard in asm_to_ref_alignment encodes both sample name and haplotype (1/2)
        # ex. HG00171_1 -> HG00171
        combined_assembly=lambda wc: os.path.join(
            config["ident_cen_ctgs"]["comb_assemblies_dir"],
            f"{wc.sm.split('_')[0]}.vrk-ps-sseq.asm-comb-dedup.fa.gz",
        ),
    output:
        fwd_cen_regions=os.path.join(
            config["ident_cen_ctgs"]["output_dir"], "{sm}_centromeric_regions.fwd.bed"
        ),
        rev_cen_regions=os.path.join(
            config["ident_cen_ctgs"]["output_dir"], "{sm}_centromeric_regions.rev.bed"
        ),
        fwd_cen_seq=os.path.join(
            config["ident_cen_ctgs"]["output_dir"], "{sm}_centromeric_regions.fwd.fa"
        ),
        rev_cen_seq=os.path.join(
            config["ident_cen_ctgs"]["output_dir"], "{sm}_centromeric_regions.rev.fa"
        ),
    log:
        "logs/extract_fwd_rev_regions_{sm}.log",
    conda:
        "../env/tools.yaml"
    shell:
        # Changed to $4 position for strand
        """
        awk -v OFS="\\t" '{{if ($6=="+") print}}' {input.all_regions} > {output.fwd_cen_regions}
        awk -v OFS="\\t" '{{if ($6=="-") print}}' {input.all_regions} > {output.rev_cen_regions}
        seqtk subseq {input.combined_assembly} {output.fwd_cen_regions} > {output.fwd_cen_seq}
        seqtk subseq {input.combined_assembly} {output.rev_cen_regions} | seqtk seq -r - > {output.rev_cen_seq}
        """


# 1. name: haplotype1-0000027
# 2. start: 96023560
# 3. end: 101450776
# 4. length: 158780978
# 5. chr: chr7
# 6. chr_coords: 58924390-64604808
# 7. orientation: +
# 8. start_end_diff: 5427216
# 9-12. results/cens/HG00171    1       centromeric     regions.fwd.bed
# haplotype1-0000027    results/cens/HG00171_chr7_haplotype1-0000027
rule create_fwd_ctg_name_legend:
    input:
        regions=rules.extract_fwd_rev_regions.output.fwd_cen_regions,
    output:
        os.path.join(config["ident_cen_ctgs"]["output_dir"], "{sm}.legend.fwd.txt"),
    log:
        "logs/create_fwd_ctg_name_legend_{sm}.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ awk -v OFS="\\t" '{{print $0, FILENAME}}' {input.regions} | \
        sed -e 's/_/\\t/g' -e 's/:/\\t/g' | \
        awk -v OFS="\\t" '{{print $1, $9"_"$5"_"$1}}' | \
        sort -k2,2;}} > {output} 2> {log}
        """


use rule create_fwd_ctg_name_legend as create_rev_ctg_name_legend with:
    input:
        regions=rules.extract_fwd_rev_regions.output.rev_cen_regions,
    output:
        os.path.join(config["ident_cen_ctgs"]["output_dir"], "{sm}.legend.rev.txt"),
    log:
        "logs/create_rev_ctg_name_legend_{sm}.log",


# :Before:
# >haplotype1-0000001:56723941-58474753
# :After:
# >
# haplotype1-0000001
# 56723941-58474753
rule split_fwd_cens_assembly_fasta:
    input:
        rules.extract_fwd_rev_regions.output.fwd_cen_seq,
    output:
        os.path.join(config["ident_cen_ctgs"]["output_dir"], "{sm}.fwd.txt"),
    conda:
        "../env/tools.yaml"
    log:
        "logs/split_fwd_cens_assembly_fasta_{sm}.log",
    shell:
        """
        sed -e 's/>/>\\n/g' -e 's/:/\\n/g' {input} > {output} 2> {log}
        """


use rule split_fwd_cens_assembly_fasta as split_rev_cens_assembly_fasta with:
    input:
        rules.extract_fwd_rev_regions.output.rev_cen_seq,
    output:
        os.path.join(config["ident_cen_ctgs"]["output_dir"], "{sm}.rev.txt"),
    log:
        "logs/split_rev_cens_assembly_fasta_{sm}.log",


# Desired: >haplotype1-0000001:chr#:6723941-58474753
rule rename_cens_fwd_ctgs:
    input:
        legend=rules.create_fwd_ctg_name_legend.output,
        seq=rules.split_fwd_cens_assembly_fasta.output,
    output:
        os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "{sm}_centromeric_regions.renamed.fwd.fa",
        ),
    conda:
        "../env/tools.yaml"
    log:
        "logs/rename_cens_fwd_ctgs_{sm}.log",
    shell:
        # Construct associative array from input legend. (line 2)
        """
        {{ awk 'BEGIN{{FS=OFS="\\t"}} \
        NR==FNR {{legend[$1]=$2; next}} \
        {{print ($1 in legend ? legend[$1] : $1)}}' {input.legend} {input.seq} | \
        awk '{{printf "%s%s", (/>/ ? ors : OFS), $0; ors=ORS}} END{{print ":"}}' | \
        sed -e 's/> />/g' -e 's/\([0-9]\) \([0-9]\)/\\1:\\2/g' | \
        tr " " "\\n";}} > {output} 2> {log}
        """


use rule rename_cens_fwd_ctgs as rename_cens_rev_ctgs with:
    input:
        legend=rules.create_rev_ctg_name_legend.output,
        seq=rules.split_rev_cens_assembly_fasta.output,
    output:
        os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "{sm}_centromeric_regions.renamed.rev.fa",
        ),
    log:
        "logs/rename_cens_rev_ctgs_{sm}.log",


rule index_renamed_cens_ctgs:
    input:
        fwd=rules.rename_cens_fwd_ctgs.output,
        rev=rules.rename_cens_rev_ctgs.output,
    output:
        fwd=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "{sm}_centromeric_regions.renamed.fwd.fai",
        ),
        rev=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "{sm}_centromeric_regions.renamed.rev.fai",
        ),
    conda:
        "../env/tools.yaml"
    log:
        "logs/index_renamed_cens_ctgs_{sm}.log",
    shell:
        """
        samtools faidx {input} &> {log}
        """


rule ident_cen_ctgs_all:
    input:
        expand(rules.intersect_cen_regions.output, sm=SAMPLES_DF.index),
        expand(rules.collapse_cens_contigs.output, sm=SAMPLES_DF.index),
        expand(
            rules.collapse_cens_contigs_only_t2t_cens.output,
            sm=SAMPLES_DF.index,
        ),
        expand(rules.intersect_filter_both_cens_ctgs.output, sm=SAMPLES_DF.index),
        expand(rules.collapse_intersected_filtered_cen_ctgs.output, sm=SAMPLES_DF.index),
        expand(rules.reintersect_sort_uniq_cens_ctgs.output, sm=SAMPLES_DF.index),
        expand(rules.extract_fwd_rev_regions.output, sm=SAMPLES_DF.index),
        # Renaming section
        expand(rules.create_fwd_ctg_name_legend.output, sm=SAMPLES_DF.index),
        expand(rules.create_rev_ctg_name_legend.output, sm=SAMPLES_DF.index),
        expand(rules.split_fwd_cens_assembly_fasta.output, sm=SAMPLES_DF.index),
        expand(rules.split_rev_cens_assembly_fasta.output, sm=SAMPLES_DF.index),
        expand(rules.rename_cens_fwd_ctgs.output, sm=SAMPLES_DF.index),
        expand(rules.rename_cens_rev_ctgs.output, sm=SAMPLES_DF.index),
        expand(rules.index_renamed_cens_ctgs.output, sm=SAMPLES_DF.index),
