
include: "common.smk"


rule filter_annotations_hor_stv:
    input:
        all_cdr_bed=lambda wc: (
            os.path.join(
                config["cdr_finder"]["output_dir"],
                "bed",
                "all_cdrs.bed",
            )
            if config.get("cdr_finder")
            else []
        ),
    output:
        cdr_bed=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "bed",
            "{chr}",
            "cdr.bed",
        ),
    params:
        cdr_output=bool(config.get("cdr_finder", False)),
    log:
        "logs/plot_hor_stv/filter_annotations_{chr}_hor_stv.log",
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


use rule modify_cenplot_tracks as modify_hor_stv_cenplot_tracks with:
    input:
        plot_layout="workflow/scripts/cenplot_hor_stv_plot.toml",
        infiles=[
            os.path.join(
                config["plot_hor_stv"]["output_dir"], "bed", "{chr}", "stv_all.bed"
            )
        ],
    output:
        plot_layout=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "plots",
            "{typ}_cens_{chr}.yaml",
        ),
    params:
        typ="{typ}",


use rule plot_multiple_cen as plot_hor_stv with:
    input:
        bed_files=[
            os.path.join(
                config["plot_hor_stv"]["output_dir"],
                "bed",
                    "{chr}",
                    "sat_annot.bed",
                ),
                os.path.join(
                config["plot_hor_stv"]["output_dir"], "bed", "{chr}", "stv_all.bed"
            ),
            rules.filter_annotations_hor_stv.output.cdr_bed,
        ],
        script="workflow/scripts/plot_multiple_cen.py",
        plot_layout=expand(
            rules.modify_hor_stv_cenplot_tracks.output, chr="{chr}", typ="all"
        ),
        # hor_stv_colors=config["plot_hor_stv"]["stv_annot_colors"],
    output:
        plots=multiext(
            os.path.join(
                config["plot_hor_stv"]["output_dir"],
                "plots",
                "all_{chr}",
            ),
            ".pdf",
            ".png",
        ),
        plot_dir=directory(
            os.path.join(
                config["plot_hor_stv"]["output_dir"],
                "plots",
                "all_{chr}",
            )
        ),
    log:
        "logs/plot_hor_stv/plot_{chr}_hor_stv_all.log",


use rule plot_multiple_cen as plot_hor_stv_complete with:
    input:
        bed_files=[
            os.path.join(
                config["plot_hor_stv"]["output_dir"],
                "bed",
                "{chr}",
                "sat_annot.bed",
            ),
            os.path.join(
                config["plot_hor_stv"]["output_dir"],
                "bed",
                "{chr}",
                "stv_complete.bed",
            ),
            rules.filter_annotations_hor_stv.output.cdr_bed,
        ],
        script="workflow/scripts/plot_multiple_cen.py",
        plot_layout=expand(
            rules.modify_hor_stv_cenplot_tracks.output, chr="{chr}", typ="complete"
        ),
    output:
        plots=multiext(
            os.path.join(
                config["plot_hor_stv"]["output_dir"],
                "plots",
                "complete_{chr}",
            ),
            ".png",
            ".pdf",
        ),
        plot_dir=directory(
            os.path.join(
                config["plot_hor_stv"]["output_dir"],
                "plots",
                "complete_{chr}",
            )
        ),
    log:
        "logs/plot_hor_stv/plot_{chr}_hor_stv_complete.log",


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
