include: "common.smk"


EXT_REF_HOR_ARR_OUTDIR = join(OUTPUT_DIR, "2-extract_ref_hor_arrays")
EXT_REF_HOR_ARR_LOGDIR = join(LOG_DIR, "2-extract_ref_hor_arrays")
EXT_REF_HOR_ARR_BMKDIR = join(BMK_DIR, "2-extract_ref_hor_arrays")


rule adjust_ref_hor_arrays:
    input:
        config["ident_cen_ctgs"]["ref_cens_regions"],
    output:
        join(
            EXT_REF_HOR_ARR_OUTDIR,
            f"{REF_NAME}.hor_arrays.{REF_CENS_EDGE_LEN}kbp.bed",
        ),
    conda:
        "../envs/tools.yaml"
    params:
        added_bases=config["extract_ref_hor_arrays"].get("added_bases", 0),
    shell:
        """
        awk -v OFS="\\t" '{{ print $1, $2-{params.added_bases}, $3+{params.added_bases} }}' {input} > {output}
        """


# Extract HOR arrays from reference.
rule extract_ref_hor_arrays:
    input:
        ref=(
            config["align_asm_to_ref"]["reference"]
            if config["align_asm_to_ref"].get("reference")
            else REF_FA
        ),
        cens_regions=rules.adjust_ref_hor_arrays.output,
    output:
        seq=join(
            EXT_REF_HOR_ARR_OUTDIR,
            f"{REF_NAME}.hor_arrays.{REF_CENS_EDGE_LEN}kbp.fa",
        ),
        idx=join(
            EXT_REF_HOR_ARR_OUTDIR,
            f"{REF_NAME}.hor_arrays.{REF_CENS_EDGE_LEN}kbp.fa.fai",
        ),
    log:
        join(EXT_REF_HOR_ARR_LOGDIR, "extract_ref_hor_arrays.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        seqtk subseq {input.ref} {input.cens_regions} > {output.seq} 2> {log}
        samtools faidx {output.seq} &> {log}
        """


rule extract_ref_hor_arrays_all:
    input:
        rules.extract_ref_hor_arrays.output,
    default_target: True
