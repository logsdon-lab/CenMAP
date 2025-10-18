# TODO: Include sample for reference at some point.


rule create_rm_bed_by_sm:
    input:
        rm_out=rules.format_repeatmasker_output.output,
    output:
        rm_bed=join(
            RM_OUTDIR,
            "bed",
            "sm_{sm}_og.bed",
        ),
    params:
        chr_rgx="",
        color_mapping=config["repeatmasker"]["repeat_colors"],
        script=workflow.source_path("../scripts/create_rm_bed.py"),
        to_abs="--to_abs",
    log:
        join(RM_LOGDIR, "create_og_rm_bed_{sm}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        shell_create_rm_bed


rule plot_og_rm_bed_by_sm:
    input:
        rm=[rules.create_rm_bed_by_sm.output],
    output:
        plots=multiext(
            join(
                RM_OUTDIR,
                "plots",
                "sm_{sm}_cens_og",
            ),
            ".pdf",
            ".png",
        ),
    wildcard_constraints:
        sm="|".join(SAMPLE_NAMES),
    params:
        script=workflow.source_path("../scripts/plot_multiple_cen.py"),
        plot_layout=workflow.source_path("../scripts/cenplot_repeatmasker_plot.yaml"),
        output_prefix=lambda wc, output: os.path.splitext(output.plots[0])[0],
        json_file_str=lambda wc, input: json.dumps(dict(input)),
        options="",
        omit_if_empty="",
        ref_ax_idx="--ref_ax_idx 0",
    log:
        join(RM_LOGDIR, "plot_og_rm_bed_{sm}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        shell_plot_multiple_cen
