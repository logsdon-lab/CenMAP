
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
        all_stv_bed=rules.aggregate_format_all_stv_row.output,
        complete_stv_bed=rules.filter_complete_correct_stv_row.output,
        all_rm_sat_bed=rules.create_rm_satellite_annotations.output,
    output:
        all_stv_bed=join(
            PLT_HOR_STV_OUTDIR,
            "bed",
            "{chr}",
            "stv_all.bed",
        ),
        complete_stv_bed=join(
            PLT_HOR_STV_OUTDIR,
            "bed",
            "{chr}",
            "stv_complete.bed",
        ),
        rm_bed=join(
            PLT_HOR_STV_OUTDIR,
            "bed",
            "{chr}",
            "sat_annot.bed",
        ),
        cdr_bed=join(
            PLT_HOR_STV_OUTDIR,
            "bed",
            "{chr}",
            "cdr.bed",
        ),
    params:
        rgx_chr=lambda wc: f"{wc.chr}[_:-]",
        cdr_output=bool(config.get("cdr_finder", False)),
    log:
        join(PLT_HOR_STV_LOGDIR, "filter_annotations_{chr}_hor_stv.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        if [ ! "{params.cdr_output}" == "False" ]; then
            {{ grep '{params.rgx_chr}' {input.all_cdr_bed} || true ;}} > {output.cdr_bed}
        fi
        {{ grep '{params.rgx_chr}' {input.all_stv_bed} || true ;}} > {output.all_stv_bed}
        {{ grep '{params.rgx_chr}' {input.complete_stv_bed} || true ;}} > {output.complete_stv_bed}
        {{ grep '{params.rgx_chr}' {input.all_rm_sat_bed} || true ;}} > {output.rm_bed}
        touch {output}
        """


rule plot_hor_stv:
    input:
        **(
            {"cdr": rules.filter_annotations_hor_stv.output.cdr_bed}
            if config.get("cdr_finder")
            else {}
        ),
        stv=lambda wc: (
            rules.filter_annotations_hor_stv.output.complete_stv_bed
            if wc.typ == "complete"
            else rules.filter_annotations_hor_stv.output.all_stv_bed
        ),
        rm=rules.filter_annotations_hor_stv.output.rm_bed,
    output:
        plots=multiext(
            join(
                PLT_HOR_STV_OUTDIR,
                "plots",
                "{typ}_{chr}",
            ),
            ".png",
            ".pdf",
        ),
    params:
        script=workflow.source_path("../scripts/plot_multiple_cen.py"),
        plot_layout=(
            workflow.source_path("../scripts/cenplot_hor_stv_plot.yaml")
            if config["humas_annot"]["mode"] != "sf"
            else workflow.source_path("../scripts/cenplot_sf_plot.yaml")
        ),
        output_prefix=lambda wc, output: os.path.splitext(output.plots[0])[0],
        json_file_str=lambda wc, input: json.dumps(dict(input)),
        options=lambda wc: f"--options '{json.dumps({"hor":{"color_map_file": os.path.abspath(config["plot_hor_stv"]["stv_annot_colors"])}})}'",
    log:
        join(PLT_HOR_STV_LOGDIR, "plot_{chr}_hor_stv_{typ}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        shell_plot_multiple_cen


rule plot_hor_stv_all:
    input:
        expand(rules.plot_hor_stv.output, chr=CHROMOSOMES, typ=["complete", "all"]),
    default_target: True
