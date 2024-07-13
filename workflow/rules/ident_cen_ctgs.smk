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
        "logs/ident_cen_ctgs/format_hor_ref_aln_cen_contigs_{sm}.log",
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
        "logs/ident_cen_ctgs/intersect_ref_cen_pqarm_{sm}.log",
    shell:
        """
        bedtools intersect -loj -a {input.aln_cens_bed} -b  {input.ref_monomeric_bed} > {output.qarms_cen_regions} 2> {log}
        """


# TODO: Don't rename sequence.
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
        resolved_cen_regions_rc=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            f"{REF_NAME}_cens",
            "{sm}_mapped_cens_rc.bed",
        ),
    params:
        thr_ctg_len=0,
    conda:
        "../env/py.yaml"
    log:
        "logs/ident_cen_ctgs/map_collapse_cens_{sm}.log",
    shell:
        """
        python {input.script} -i {input.regions} -t {params.thr_ctg_len} > {output.resolved_cen_regions} 2> {log}
        awk -v OFS="\\t" '{{
            if ($6 == "true") {{
                $4="rc-"$4
            }};print
        }}' {output.resolved_cen_regions} > {output.resolved_cen_regions_rc} 2>> {log}
        """


# Make renamed copy of assembly here with mapped chr.
RENAME_ASM_CFG = {
    "bed_input_regions": rules.map_collapse_cens.output.resolved_cen_regions_rc,
    "fa_assembly": os.path.join(
        config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fasta"
    ),
    "output_dir": os.path.join(config["concat_asm"]["output_dir"], "{sm}"),
    "samples": SAMPLE_NAMES,
    "logs_dir": "logs/ident_cen_ctgs/rename_cens",
}


module rename_asm:
    snakefile:
        "rename_ctgs.smk"
    config:
        RENAME_ASM_CFG


use rule * from rename_asm as asm_*


# TODO:
rule fix_ort_asm:
    input:
        # This is the key. Update as go along.
        bed=rules.asm_create_renamed_bed_n_legend.output.regions_renamed,
        fa=rules.asm_rename_ctgs.output,
        fai=rules.asm_index_renamed_ctgs.output,
    output:
        fa=os.path.join(
            config["concat_asm"]["output_dir"], "{sm}", "{sm}_regions.renamed.reort.fa"
        ),
        faidx=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}",
            "{sm}_regions.renamed.reort.fa.fai",
        ),
    log:
        "logs/ident_cen_ctgs/fix_ort_asm_{sm}.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        # Reverse complement sequences and then get everything else.
        cat <(seqtk subseq <(seqtk seq -r {input.fa}) <(grep "rc-" {input.bed} | cut -f 1)) \
            <(seqtk subseq {input.fa} <(grep -v -f <(grep "rc-" {input.bed} | cut -f 1) {input.fai} | cut -f 1)) \
        > {output.fa} 2> {log}
        samtools faidx {output.fa} 2> {log}
        """


use rule extract_and_index_fa as extract_cens_regions with:
    input:
        bed=rules.asm_create_renamed_bed_n_legend.output.regions_renamed,
        fa=rules.fix_ort_asm.output.fa,
    output:
        seq=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "seq",
            "{sm}_centromeric_regions.fa",
        ),
        idx=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "seq",
            "{sm}_centromeric_regions.fa.fai",
        ),
    log:
        "logs/ident_cen_ctgs/extract_regions_{sm}.log",
    conda:
        "../env/tools.yaml"


rule ident_cen_ctgs_all:
    input:
        # Rename ctgs
        rules.asm_rename_ctg_all.input,
        # Fix orientation.
        expand(rules.fix_ort_asm.output, sm=SAMPLE_NAMES),
        expand(
            rules.extract_cens_regions.output,
            sm=SAMPLE_NAMES,
        ),
