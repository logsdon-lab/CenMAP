# Identify and add sample and chromosome to centromeric regions or just add sample.


include: "common.smk"


IDENT_CEN_CTGS_OUTDIR = join(OUTPUT_DIR, "5-ident_cen_ctgs")
IDENT_CEN_CTGS_LOGDIR = join(LOG_DIR, "5-ident_cen_ctgs")
IDENT_CEN_CTGS_BMKDIR = join(BMK_DIR, "5-ident_cen_ctgs")


if CHROMOSOMES:

    include: "5.1-add_sm_chr_cen_ctgs.smk"

else:

    include: "5.1-add_sm_cen_ctgs.smk"


# Use SRF regions as putative alpha-satellite regions.
rule make_srf_putative_alr_regions:
    input:
        bed=rules.extract_filter_monomers.output.bed_mon,
        rename_key=rules.create_rename_key.output if CHROMOSOMES else [],
    output:
        join(
            IDENT_CEN_CTGS_OUTDIR,
            "bed",
            "{sm}_alr.bed",
        ),
    params:
        bp_min_length=config["ident_cen_ctgs"]["bp_min_length"],
        # Strict merging (Only 1 LINE element) as we need this to be large blocks.
        bp_merge=8000,
        format_cmd=lambda wc, input: (
            " ".join(
                    [
                        "join",
                        "-",
                        f"<(sort -k1,1 {input.rename_key})",
                        "|",
                        "awk",
                        "-v",
                        "FS=' '",
                        "-v",
                        "OFS='\\t'",
                        """'{{ if ($10 ~ "rc-chr") {{ st=$11-$3; end=$11-$2; }} else {{ st=$2; end=$3; }}; print $10, st, end, "ALR/Alpha" }}'""",
                        "|",
                        "sort",
                        "|",
                        "uniq",
                    ]
                )
                if input.rename_key
            else f"""awk -v OFS="\\t" '{{{{ $1="{wc.sm}_"$1; print $1, $2, $3, "ALR/Alpha" }}}}'"""
        ),
    log:
        join(IDENT_CEN_CTGS_LOGDIR, "make_srf_putative_alr_regions_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ srf-n-trf regions \
        -b {input.bed} \
        -d {params.bp_merge} \
        -m {params.bp_min_length} | \
        sort -k1,1 | \
        {params.format_cmd} ;}} > {output} 2> {log}
        """
