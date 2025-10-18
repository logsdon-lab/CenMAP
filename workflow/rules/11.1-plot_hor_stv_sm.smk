
rule filter_annotations_hor_stv_by_sm:
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
        all_stv_bed=expand(
            rules.aggregate_format_all_stv_row.output,
            chr=CHROMOSOMES if CHROMOSOMES else ["all"],
        ),
        complete_stv_bed=expand(
            rules.filter_complete_correct_stv_row.output,
            chr=CHROMOSOMES if CHROMOSOMES else ["all"],
        ),
        all_rm_sat_bed=rules.create_rm_satellite_annotations.output,
    output:
        all_stv_bed=join(
            PLT_HOR_STV_OUTDIR,
            "bed",
            "sm_stv_all_{sm}.bed",
        ),
        complete_stv_bed=join(
            PLT_HOR_STV_OUTDIR,
            "bed",
            "sm_stv_complete_{sm}.bed",
        ),
        rm_bed=join(
            PLT_HOR_STV_OUTDIR,
            "bed",
            "sm_sat_annot_{sm}.bed",
        ),
        cdr_bed=join(
            PLT_HOR_STV_OUTDIR,
            "bed",
            "sm_cdr_{sm}.bed",
        ),
    wildcard_constraints:
        sm="|".join(SAMPLE_NAMES),
    params:
        rgx=lambda wc: f"^{wc.sm}[_:-]",
        cdr_output=bool(config.get("cdr_finder", False)),
    log:
        join(PLT_HOR_STV_LOGDIR, "filter_annotations_{sm}_hor_stv.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        if [ ! "{params.cdr_output}" == "False" ]; then
            {{ grep -P '{params.rgx}' {input.all_cdr_bed} || true ;}} > {output.cdr_bed}
        fi
        {{ grep -P '{params.rgx}' {input.all_stv_bed} || true ;}} > {output.all_stv_bed}
        {{ grep -P '{params.rgx}' {input.complete_stv_bed} || true ;}} > {output.complete_stv_bed}
        {{ grep -P '{params.rgx}' {input.all_rm_sat_bed} || true ;}} > {output.rm_bed}
        """


rule plot_hor_stv_by_sm:
    input:
        **(
            {"cdr": rules.filter_annotations_hor_stv_by_sm.output.cdr_bed}
            if config.get("cdr_finder")
            else {}
        ),
        stv=lambda wc: (
            rules.filter_annotations_hor_stv_by_sm.output.complete_stv_bed
            if wc.typ == "complete"
            else rules.filter_annotations_hor_stv_by_sm.output.all_stv_bed
        ),
        rm=rules.filter_annotations_hor_stv_by_sm.output.rm_bed,
    output:
        plots=multiext(
            join(
                PLT_HOR_STV_OUTDIR,
                "plots",
                "sm_{typ}_{sm}",
            ),
            ".png",
            ".pdf",
        ),
    wildcard_constraints:
        sm="|".join(SAMPLE_NAMES),
    params:
        script=workflow.source_path("../scripts/plot_multiple_cen.py"),
        plot_layout=get_hor_stv_plot_layout,
        output_prefix=lambda wc, output: os.path.splitext(output.plots[0])[0],
        json_file_str=lambda wc, input: json.dumps(dict(input)),
        options=lambda wc: f"--options '{json.dumps({"hor":{"color_map_file": os.path.abspath(config["plot_hor_stv"]["stv_annot_colors"])}})}'",
        omit_if_empty=lambda wc: "--omit_if_any_empty" if wc.typ == "complete" else "",
        ref_ax_idx="--ref_ax_idx 1",
    log:
        join(PLT_HOR_STV_LOGDIR, "plot_{sm}_hor_stv_{typ}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        shell_plot_multiple_cen
