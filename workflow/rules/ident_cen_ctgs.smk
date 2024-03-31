# Identify centromeric regions
# Extract centromeric regions


# Get centromeric regions from alignments to t2t-chm13 ONLY
rule intersect_cen_regions:
    input:
        aln_bed=lambda wc: expand(
            rules.asm_ref_aln_to_bed.output, ref=[REF_NAME], sm=wc.sm
        ),
        ref_cens_bed=config["ident_cen_ctgs"]["ref_cens_500kbp_regions"],
    output:
        cen_regions=os.path.join(
            config["ident_cen_ctgs"]["output_dir"], "bed", "{sm}_cens.bed"
        ),
    conda:
        "../env/tools.yaml"
    log:
        "logs/intersect_ref_cen_{sm}.log",
    shell:
        """
        bedtools intersect -a {input.aln_bed} -b {input.ref_cens_bed} > {output} 2> {log}
        """


rule intersect_with_pq_arm:
    input:
        aln_cens_bed=rules.intersect_cen_regions.output,
        ref_monomeric_bed=config["ident_cen_ctgs"]["ref_cens_monomeric_regions"],
    output:
        pqarms_cen_regions=os.path.join(
            config["ident_cen_ctgs"]["output_dir"], "bed", "{sm}_pqarm_cens.bed"
        ),
    conda:
        "../env/tools.yaml"
    log:
        "logs/intersect_ref_cen_pqarm_{sm}.log",
    shell:
        """
        bedtools intersect -a {input.aln_cens_bed} -b {input.ref_monomeric_bed} > {output} 2> {log}
        """


# In.
# 6. query_name
# 7. query_start
# 8. query_end
# 1. reference_name
# 5. strand
# Out.
# 1. query_name
# 2. query_start
# 3. query_end
# 5. reference_name
# 6. strand
# +. sub(query_end, query_start)
rule collapse_cens_contigs:
    input:
        script="workflow/scripts/map_cens.py",
        regions=rules.intersect_with_pq_arm.output,
    output:
        regions=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            "{sm}_centromeric_contigs.bed",
        ),
    params:
        len_thr=1_000_000,
    conda:
        "../env/py.yaml"
    log:
        "logs/collapse_cent_ctgs_{sm}.log",
    shell:
        """
        python {input.script} -i {input.regions} -t {params.len_thr} > {output.regions} 2> {log}
        """


# Make renamed copy of assembly here with mapped chr.
RENAME_ASM_CFG = {
    "bed_input_regions": rules.collapse_cens_contigs.output,
    "fa_assembly": rules.concat_asm.output,
    "output_dir": os.path.join(config["concat_asm"]["output_dir"], "{sm}"),
    "samples": SAMPLE_NAMES,
    "log_dir": "logs/rename_cens",
}


module rename_asm:
    snakefile:
        "rename_ctgs.smk"
    config:
        RENAME_ASM_CFG


use rule * from rename_asm as asm_*


rule filter_cens_oriented_regions:
    input:
        all_regions=lambda wc: expand(
            rules.asm_create_renamed_bed_n_legend.output.regions_renamed,
            ort=ORIENTATION,
            sm=[wc.sm],
        ),
    output:
        regions=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            "{sm}_centromeric_regions.{ort}.bed",
        ),
    params:
        sign=lambda wc: "+" if wc.ort == "fwd" else "-",
    log:
        "logs/filter_{ort}_regions_{sm}.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        awk -v OFS="\\t" '{{if ($5=="{params.sign}") print}}' {input.all_regions} > {output.regions} 2> {log}
        """


use rule extract_and_index_fa as extract_cens_oriented_regions with:
    input:
        bed=rules.filter_cens_oriented_regions.output.regions,
        fa=rules.asm_rename_ctgs.output,
    output:
        seq=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "seq",
            "{sm}_centromeric_regions.{ort}.fa",
        ),
        idx=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "seq",
            "{sm}_centromeric_regions.{ort}.fa.fai",
        ),
    params:
        added_cmds=lambda wc: "" if wc.ort == "fwd" else "| seqtk seq -r",
    log:
        "logs/extract_{ort}_regions_{sm}.log",


rule ident_cen_ctgs_all:
    input:
        expand(rules.intersect_cen_regions.output, sm=SAMPLE_NAMES),
        expand(rules.intersect_with_pq_arm.output, sm=SAMPLE_NAMES),
        expand(rules.collapse_cens_contigs.output, sm=SAMPLE_NAMES),
        # Rename ctgs
        rules.asm_rename_ctg_all.input,
        expand(
            rules.filter_cens_oriented_regions.output,
            sm=SAMPLE_NAMES,
            ort=ORIENTATION,
        ),
        expand(
            rules.extract_cens_oriented_regions.output,
            sm=SAMPLE_NAMES,
            ort=ORIENTATION,
        ),
