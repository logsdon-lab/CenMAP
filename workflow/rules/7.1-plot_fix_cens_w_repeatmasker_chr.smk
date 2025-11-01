
rule create_fixed_rm_bed_by_chr:
    input:
        rm_out=expand(rules.fix_cens_rm_out.output, sm=SAMPLE_NAMES),
    output:
        rm_bed=join(
            FIX_RM_OUTDIR,
            "bed",
            "chr_{chr}.bed",
        ),
    params:
        chr_rgx=lambda wc: f"-c {wc.chr}[:_-]" if wc.chr != "all" else "",
        color_mapping=config["repeatmasker"]["repeat_colors"],
        script=workflow.source_path("../scripts/create_rm_bed.py"),
        to_abs="",
    log:
        join(FIX_RM_LOGDIR, "create_fixed_rm_bed_{chr}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        shell_create_rm_bed


rule plot_fixed_rm_bed_by_chr:
    input:
        rm=[
            (
                rules.create_ref_rm_bed_by_chr.output
                if config["repeatmasker"]["ref_repeatmasker_output"]
                else []
            ),
            rules.create_fixed_rm_bed_by_chr.output,
        ],
    output:
        plots=multiext(
            join(
                FIX_RM_OUTDIR,
                "plots",
                "chr_{chr}_cens",
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
        join(FIX_RM_LOGDIR, "plot_fixed_rm_bed_{chr}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        shell_plot_multiple_cen
