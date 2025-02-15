
include: "common.smk"
include: "utils.smk"
# include: "8-cdr_finder.smk"
include: "8-format_repeatmasker_sat_annot.smk"
include: "10-format_hor_stv.smk"


PLT_HOR_STV_OUTDIR = join(OUTPUT_DIR, "11-plot_hor_stv")
PLT_HOR_STV_LOGDIR = join(LOG_DIR, "11-plot_hor_stv")
PLT_HOR_STV_BMKDIR = join(BMK_DIR, "11-plot_hor_stv")


rule filter_annotations_hor_stv:
    input:
        all_cdr_bed=(
            join(
                OUTPUT_DIR,
                "8-cdr_finder",
                "bed",
                "all_cdrs.bed",
            )
            if config.get("cdr_finder")
            else []
        ),
    output:
        cdr_bed=join(
            PLT_HOR_STV_OUTDIR,
            "bed",
            "{chr}",
            "cdr.bed",
        ),
    params:
        cdr_output=bool(config.get("cdr_finder", False)),
    log:
        join(PLT_HOR_STV_LOGDIR, "filter_annotations_{chr}_hor_stv.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        if [ "{params.cdr_output}" == "False" ]; then
            touch {output.cdr_bed}
        else
            ( grep '{wildcards.chr}_' {input.all_cdr_bed} || true ) > {output.cdr_bed}
        fi
        """


# Require both to exist before beginning plotting.
rule modify_hor_stv_cenplot_tracks:
    input:
        plot_layout=workflow.source_path("../scripts/cenplot_hor_stv_plot.toml"),
        infiles=[
            rules.aggregate_format_all_stv_row.output,
            rules.filter_complete_correct_stv_row.output,
            rules.split_rm_satellite_annotations.output,
        ],
    output:
        plot_layout=join(
            PLT_HOR_STV_OUTDIR,
            "plots",
            "{typ}_cens_{chr}.yaml",
        ),
    params:
        typ="{typ}",
    run:
        format_toml_path(
            input_plot_layout=input.plot_layout,
            output_plot_layout=output.plot_layout,
            indir=os.path.abspath(os.path.dirname(str(input.infiles[0]))),
            **dict(params.items()),
        )


rule plot_hor_stv:
    input:
        bed_files=[
            rules.split_rm_satellite_annotations.output,
            rules.aggregate_format_all_stv_row.output,
            rules.filter_annotations_hor_stv.output.cdr_bed,
        ],
        script=workflow.source_path("../scripts/plot_multiple_cen.py"),
        plot_layout=expand(
            rules.modify_hor_stv_cenplot_tracks.output, chr="{chr}", typ="all"
        ),
        # hor_stv_colors=config["plot_hor_stv"]["stv_annot_colors"],
    output:
        plots=multiext(
            join(
                PLT_HOR_STV_OUTDIR,
                "plots",
                "all_{chr}",
            ),
            ".pdf",
            ".png",
        ),
        plot_dir=directory(
            join(
                PLT_HOR_STV_OUTDIR,
                "plots",
                "all_{chr}",
            )
        ),
    log:
        join(PLT_HOR_STV_LOGDIR, "plot_{chr}_hor_stv_all.log"),
    threads: 4
    conda:
        "../envs/py.yaml"
    shell:
        shell_plot_multiple_cen


rule plot_hor_stv_complete:
    input:
        bed_files=[
            rules.split_rm_satellite_annotations.output,
            rules.filter_complete_correct_stv_row.output,
            rules.filter_annotations_hor_stv.output.cdr_bed,
        ],
        script=workflow.source_path("../scripts/plot_multiple_cen.py"),
        plot_layout=expand(
            rules.modify_hor_stv_cenplot_tracks.output, chr="{chr}", typ="complete"
        ),
    output:
        plots=multiext(
            join(
                PLT_HOR_STV_OUTDIR,
                "plots",
                "complete_{chr}",
            ),
            ".png",
            ".pdf",
        ),
        plot_dir=directory(
            join(
                PLT_HOR_STV_OUTDIR,
                "plots",
                "complete_{chr}",
            )
        ),
    log:
        join(PLT_HOR_STV_LOGDIR, "plot_{chr}_hor_stv_complete.log"),
    threads: 4
    conda:
        "../envs/py.yaml"
    shell:
        shell_plot_multiple_cen


rule plot_hor_stv_all:
    input:
        expand(
            rules.plot_hor_stv.output,
            chr=CHROMOSOMES,
        ),
        expand(
            rules.plot_hor_stv_complete.output,
            chr=CHROMOSOMES,
        ),
    default_target: True
