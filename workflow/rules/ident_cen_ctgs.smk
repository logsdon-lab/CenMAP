# Identify centromeric regions
# Extract centromeric regions


include: "utils.smk"
include: "common.smk"


rule format_hor_ref_aln_cen_contigs:
    input:
        aln_bed=os.path.join("results", f"{REF_NAME}_cens", "bed", "{sm}.bed"),
    output:
        cen_regions=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            f"{REF_NAME}_cens",
            "{sm}_cens.bed",
        ),
    conda:
        "../env/tools.yaml"
    log:
        "logs/format_hor_ref_aln_cen_contigs_{sm}.log",
    shell:
        # 1. chr2:91797392-95576642
        # 2. 3054999
        # 3. 3779251
        # 4. 3779251
        # 5. -
        # 6. h1tg000001l#1-110442987
        # 7. 15032098
        # 8. 15756783
        # 9. 110442987
        # 10. 99.97791
        # 11. 99.97032
        # 12. 99.89915
        # 13. 724023
        # 14. 160
        # 15. 27
        # 16. 28
        # 17. 69
        # 18. 502
        """
        awk -v OFS="\\t" 'NR>1 {{
            # Find starts/ends in contig name.
            match($1, ":(.+)-", starts);
            # Remove coords from ctg name
            gsub(":.*-.*", "", $1)
            # Print columns.
            print $1, $2+starts[1], $3+starts[1], $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18
        }}' {input} > {output} 2> {log}
        """


# Add monomeric pq arms and intersect with alignments.
rule intersect_with_pq_arm:
    input:
        aln_cens_bed=rules.format_hor_ref_aln_cen_contigs.output,
        ref_monomeric_bed=config["ident_cen_ctgs"]["ref_cens_monomeric_regions"],
    output:
        qarms_cen_regions=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            f"{REF_NAME}_cens",
            "{sm}_pqarm_cens.bed",
        ),
    conda:
        "../env/tools.yaml"
    log:
        "logs/intersect_ref_cen_pqarm_{sm}.log",
    shell:
        """
        bedtools intersect -loj -a {input.aln_cens_bed} -b  {input.ref_monomeric_bed} > {output.qarms_cen_regions} 2> {log}
        """


rule map_collapse_cens:
    input:
        script="workflow/scripts/map_cens.py",
        regions=rules.intersect_with_pq_arm.output,
    output:
        resolved_cen_regions=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            f"{REF_NAME}_cens",
            "{sm}_mapped_cens.bed",
        ),
    params:
        thr_ctg_len=1_000_000,
    conda:
        "../env/py.yaml"
    log:
        "logs/map_collapse_cens_{sm}.log",
    shell:
        """
        python {input.script} -i {input.regions} -t {params.thr_ctg_len} > {output} 2> {log}
        """


# Make renamed copy of assembly here with mapped chr.
RENAME_ASM_CFG = {
    "bed_input_regions": rules.map_collapse_cens.output,
    "fa_assembly": os.path.join(
        config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fasta"
    ),
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
