
include: "utils.smk"
include: "common.smk"
include: "2-concat_asm.smk"
include: "3-srf.smk"


IDENT_CEN_CTGS_OUTDIR = join(OUTPUT_DIR, "5-ident_cen_ctgs")
IDENT_CEN_CTGS_LOGDIR = join(LOG_DIR, "5-ident_cen_ctgs")
IDENT_CEN_CTGS_BMKDIR = join(BMK_DIR, "5-ident_cen_ctgs")


rule create_rename_key:
    input:
        idx=rules.concat_asm.output.idx,
    output:
        # (old_name, new_name, length)
        rename_key=join(
            IDENT_CEN_CTGS_OUTDIR,
            "bed",
            "{sm}_rename_key.tsv",
        ),
    conda:
        "../envs/tools.yaml"
    log:
        join(IDENT_CEN_CTGS_LOGDIR, "create_rename_key_{sm}.log"),
    shell:
        """
        awk -v OFS="\\t" '{{ print $1, "{wildcards.sm}_"$1, $2}}' {input} > {output.rename_key} 2> {log}
        """


rule create_final_asm:
    input:
        fa=rules.concat_asm.output.fa,
        idx=rules.concat_asm.output.idx,
        rename_key=rules.create_rename_key.output,
    output:
        fa=join(
            CONCAT_ASM_OUTDIR,
            "{sm}-asm-renamed.fa",
        ),
        idx=join(
            CONCAT_ASM_OUTDIR,
            "{sm}-asm-renamed.fa.fai",
        ),
    params:
        pattern=r"'^(\S+)\s*'",
        replacement=lambda wc: "'{kv}'",
    conda:
        "../envs/tools.yaml"
    log:
        join(IDENT_CEN_CTGS_BMKDIR, "create_{sm}_final_asm.log"),
    shell:
        """
        seqkit replace -p {params.pattern} -r {params.replacement} \
        -k {input.rename_key} \
        {input.fa} \
        --keep-key > {output.fa} 2> {log}
        samtools faidx {output.fa} 2>> {log}
        """


# Reorient satellite bed.
# Don't wait for rename assembly as can do in parallel.
rule create_final_satellite_bed:
    input:
        bed=rules.merge_slop_region_bed.output,
    output:
        bed=join(
            IDENT_CEN_CTGS_OUTDIR,
            "bed",
            "{sm}.bed",
        ),
    conda:
        "../envs/tools.yaml"
    log:
        join(IDENT_CEN_CTGS_LOGDIR, "create_final_satellite_bed_{sm}.log"),
    shell:
        """
        awk -v OFS="\\t" '{{ $1="{wildcards.sm}_"$1; print }}' {input.bed} > {output.bed} 2> {log}
        """


# Extract complete centromeric contigs.
rule extract_cens_regions:
    input:
        bed=rules.create_final_satellite_bed.output,
        fa=rules.create_final_asm.output.fa,
    output:
        seq=temp(
            join(
                IDENT_CEN_CTGS_OUTDIR,
                "seq",
                "{sm}_satellite_regions.fa",
            )
        ),
        idx=temp(
            join(
                IDENT_CEN_CTGS_OUTDIR,
                "seq",
                "{sm}_satellite_regions.fa.fai",
            )
        ),
    params:
        bed=lambda wc, input: input.bed,
        added_cmds="",
    log:
        join(IDENT_CEN_CTGS_LOGDIR, "extract_regions_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        shell_extract_and_index_fa


rule ident_cen_ctgs_all:
    input:
        expand(
            rules.create_rename_key.output,
            sm=SAMPLE_NAMES,
        ),
        expand(
            rules.create_final_asm.output,
            sm=SAMPLE_NAMES,
        ),
    default_target: True
