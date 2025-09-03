include: "common.smk"
include: "7-fix_cens_w_repeatmasker.smk"


FMT_RM_SAT_OUTDIR = join(OUTPUT_DIR, "8-format_repeatmasker_sat_annot")
FMT_RM_SAT_LOGDIR = join(LOG_DIR, "8-format_repeatmasker_sat_annot")
FMT_RM_SAT_BMKDIR = join(BMK_DIR, "8-format_repeatmasker_sat_annot")


rule create_rm_satellite_annotations:
    input:
        rm_output=[
            expand(
                rules.fix_cens_rm_out.output,
                sm=SAMPLE_NAMES,
            ),
            config["repeatmasker"]["ref_repeatmasker_output"],
        ],
        rm_sat_patterns=config.get("plot_hor_stv", {}).get("sat_annot_colors", []),
    output:
        join(
            FMT_RM_SAT_OUTDIR,
            "bed",
            "all_cens.annotation.bed",
        ),
    params:
        sat_annot_colors=lambda wc, input: (
            f"--patterns {input.rm_sat_patterns}" if input.rm_sat_patterns else ""
        ),
        script=workflow.source_path("../scripts/format_rm_sat_annot.py"),
    log:
        join(FMT_RM_SAT_LOGDIR, "create_rm_satellite_annotations.log"),
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {params.script} \
        -i <(awk -v OFS="\\t" '{{$1=$1; print}}' {input.rm_output} | cut -f 1-15) \
        {params.sat_annot_colors} \
        --add_ct > {output} 2> {log}
        """


rule format_repeatmasker_sat_all:
    input:
        rules.create_rm_satellite_annotations.output,
    default_target: True
