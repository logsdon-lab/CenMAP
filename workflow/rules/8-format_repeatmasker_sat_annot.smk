include: "common.smk"
include: "7-finalize_cens.smk"


FMT_RM_SAT_OUTDIR = join(OUTPUT_DIR, "8-format_repeatmasker_sat_annot")
FMT_RM_SAT_LOGDIR = join(LOG_DIR, "8-format_repeatmasker_sat_annot")
FMT_RM_SAT_BMKDIR = join(BMK_DIR, "8-format_repeatmasker_sat_annot")


rule reformat_srf_putative_alr_to_rm:
    input:
        bed=rules.make_srf_putative_alr_regions.output,
        rename_key=rules.make_complete_cens_bed.output.rename_key,
    output:
        pipe(
            join(
                FMT_RM_SAT_OUTDIR,
                "bed",
                "{sm}.rm.annotation.srf.out",
            )
        ),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        join {input.bed} \
            <(awk -v OFS="\\t" '{{ split($1, chrom_coords, ":"); print chrom_coords[1], $2 }}' {input.rename_key}) | \
        awk -v FS=" " -v OFS="\\t" '{{
            new_name=$5
            split(new_name, chrom_coords, ":")
            split(chrom_coords[2], coords, "-")
            if ($2 >= coords[1] && $3 <= coords[2]) {{
                # [0,4,5,6,8,9,10] - columns for RM
                print 0, ".", ".", ".", new_name, $2, $3, ".", "+", $4, "."
            }}
        }}' > {output}
        """


rule create_rm_satellite_annotations:
    input:
        rm_output=(
            [
                expand(
                    rules.fix_cens_rm_out.output,
                    sm=SAMPLE_NAMES,
                ),
                config["repeatmasker"]["ref_repeatmasker_output"],
            ]
            if RUN_REPEATMASKER
            else expand(
                rules.reformat_srf_putative_alr_to_rm.output,
                sm=SAMPLE_NAMES,
            )
        ),
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
