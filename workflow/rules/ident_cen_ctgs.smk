# Identify centromeric regions
# Extract centromeric regions

OUTPUT_DIR_T2T_REF_REGIONS = f"results/{REF_NAME}/bed"


# Get centromeric regions from alignments to t2t-chm13 ONLY
rule intersect_cen_regions:
    input:
        left=lambda wc: expand(
            rules.asm_ref_aln_to_bed.output, ref=[REF_NAME], sm=wc.sm
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
        script="workflow/scripts/bedminmax.py",
        regions=rules.intersect_cen_regions.output,
    output:
        refmt_regions=temp(
            os.path.join(
                OUTPUT_DIR_T2T_REF_REGIONS,
                "len_calc_{sm}_centromeric_contigs.bed",
            )
        ),
        regions=os.path.join(OUTPUT_DIR_T2T_REF_REGIONS, "{sm}_centromeric_contigs.bed"),
    params:
        len_thr=1_000_000
    conda:
        "../env/py.yaml"
    log:
        "logs/collapse_cent_ctgs_{sm}.log",
    shell:
        """
        awk -v OFS="\\t" '{{ print $6, $7, $8, $1, $5 }}' \
        {input.regions} > {output.refmt_regions} 2> {log}

        # Calculate length of ref region.
        {{ python {input.script} -i {output.refmt_regions} | \
        awk -v OFS="\\t" '{{
            contig_len=$3-$2
            if (contig_len > {params.len_thr} ) {{ print $1, $2, $3, $4, $5, contig_len }}
        }}' |  \
        sort -k1,1;}} > {output.regions} 2> {log}
        """


# TODO: Filter 1_000_000 thr.
rule filter_cens_oriented_regions:
    input:
        all_regions=lambda wc: expand(
            rules.collapse_cens_contigs.output.regions,
            ort=ORIENTATION,
            sm=[wc.sm],
        ),
    output:
        regions=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
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
        awk -v OFS="\\t" '{{if ($6=="{params.sign}") print}}' {input.all_regions} > {output.regions} 2> {log}
        """


use rule extract_and_index_fa as extract_cens_oriented_regions with:
    input:
        bed=rules.filter_cens_oriented_regions.output,
        fa=rules.concat_asm.output,
    output:
        seq=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "{sm}_centromeric_regions.{ort}.fa",
        ),
        idx=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "{sm}_centromeric_regions.{ort}.fa.fai",
        ),
    params:
        added_cmds=lambda wc: "" if wc.ort == "fwd" else "| seqtk seq -r",
    log:
        "logs/extract_{ort}_regions_{sm}.log",


RENAME_CTGS_CFG = {
    "bed_input_regions": rules.filter_cens_oriented_regions.output.regions,
    "fa_assembly": rules.extract_cens_oriented_regions.output.seq,
    "output_dir": config["ident_cen_ctgs"]["output_dir"],
    "samples": SAMPLE_NAMES,
    "log_dir": "logs/rename_cens",
    "sed_cmd": "sed -e 's/> />/g' -e 's/\([0-9]\) \([0-9]\)/\\1:\\2/g' | tr \" \" \"\\n\"",
}


module rename_cens_ctgs:
    snakefile:
        "rename_ctgs.smk"
    config:
        RENAME_CTGS_CFG


use rule * from rename_cens_ctgs as cens_*


rule ident_cen_ctgs_all:
    input:
        expand(rules.intersect_cen_regions.output, sm=SAMPLE_NAMES),
        expand(rules.collapse_cens_contigs.output, sm=SAMPLE_NAMES),
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
        # Rename cens ctgs
        rules.cens_rename_ctg_all.input,
