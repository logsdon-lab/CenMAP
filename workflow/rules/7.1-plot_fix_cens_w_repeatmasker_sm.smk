

rule create_fixed_rm_bed_by_sm:
    input:
        rm_out=rules.fix_cens_rm_out.output,
    output:
        rm_bed=join(
            FIX_RM_OUTDIR,
            "bed",
            "sm_{sm}.bed",
        ),
    params:
        chr_rgx="",
        color_mapping=config["repeatmasker"]["repeat_colors"],
        script=workflow.source_path("../scripts/create_rm_bed.py"),
        to_abs="",
    log:
        join(FIX_RM_LOGDIR, "create_fixed_rm_bed_{sm}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        shell_create_rm_bed


rule plot_fixed_rm_bed_by_sm:
    input:
        rm=[
            rules.create_fixed_rm_bed_by_sm.output,
        ],
    output:
        plots=multiext(
            join(
                FIX_RM_OUTDIR,
                "plots",
                "sm_{sm}_cens",
            ),
            ".pdf",
            ".png",
        ),
    params:
        script=workflow.source_path("../scripts/plot_multiple_cen.py"),
        plot_layout=workflow.source_path("../scripts/cenplot_repeatmasker_plot.yaml"),
        output_prefix=lambda wc, output: os.path.splitext(output.plots[0])[0],
        json_file_str=lambda wc, input: json.dumps(dict(input)),
        options="",
        omit_if_empty="",
        ref_ax_idx="--ref_ax_idx 0",
    log:
        join(FIX_RM_LOGDIR, "plot_fixed_rm_bed_{sm}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        shell_plot_multiple_cen
