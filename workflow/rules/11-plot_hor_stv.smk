
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
        rgx_chr=lambda wc: f"{wc.chr}[_:-]" if wc.chr != "all" else "",
        cdr_output=bool(config.get("cdr_finder", False)),
    log:
        join(PLT_HOR_STV_LOGDIR, "filter_annotations_{chr}_hor_stv.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        if [ -z "{params.rgx_chr}" ]; then
            if [ ! "{params.cdr_output}" == "False" ]; then
                cp {input.all_cdr_bed} {output.cdr_bed}
            fi
            cp {input.all_stv_bed} {output.all_stv_bed}
            cp {input.complete_stv_bed} {output.complete_stv_bed}
            cp {input.all_rm_sat_bed} {output.rm_bed}
        else
            if [ ! "{params.cdr_output}" == "False" ]; then
                {{ grep '{params.rgx_chr}' {input.all_cdr_bed} || true ;}} > {output.cdr_bed}
            fi
            {{ grep '{params.rgx_chr}' {input.all_stv_bed} || true ;}} > {output.all_stv_bed}
            {{ grep '{params.rgx_chr}' {input.complete_stv_bed} || true ;}} > {output.complete_stv_bed}
            {{ grep '{params.rgx_chr}' {input.all_rm_sat_bed} || true ;}} > {output.rm_bed}
        fi
        touch {output}
        """


def get_hor_stv_plot_layout(wc) -> str:
    if config["humas_annot"]["mode"] == "sf":
        return workflow.source_path("../scripts/cenplot_sf_plot.yaml")
    elif config["humas_annot"]["mode"] == "srf-n-trf":
        return workflow.source_path("../scripts/cenplot_srf-n-trf_plot.yaml")
    else:
        return workflow.source_path("../scripts/cenplot_hor_stv_plot.yaml")


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
        plot_layout=get_hor_stv_plot_layout,
        output_prefix=lambda wc, output: os.path.splitext(output.plots[0])[0],
        json_file_str=lambda wc, input: json.dumps(dict(input)),
        options=lambda wc: f"--options '{json.dumps({"hor":{"color_map_file": os.path.abspath(config["plot_hor_stv"]["stv_annot_colors"])}})}'",
        omit_if_empty=lambda wc: "--omit_if_any_empty" if wc.typ == "complete" else "",
        ref_ax_idx="--ref_ax_idx 1",
    log:
        join(PLT_HOR_STV_LOGDIR, "plot_{chr}_hor_stv_{typ}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        shell_plot_multiple_cen


rule plot_hor_stv_all:
    input:
        expand(
            rules.plot_hor_stv.output,
            chr=CHROMOSOMES if CHROMOSOMES else ["all"],
            typ=["complete", "all"],
        ),
    default_target: True
