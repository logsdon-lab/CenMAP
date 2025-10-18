
include: "common.smk"
include: "utils.smk"
# include: "8-cdr_finder.smk"
include: "8-format_repeatmasker_sat_annot.smk"
include: "10-format_hor_stv.smk"


PLT_HOR_STV_OUTDIR = join(OUTPUT_DIR, "11-plot_hor_stv")
PLT_HOR_STV_LOGDIR = join(LOG_DIR, "11-plot_hor_stv")
PLT_HOR_STV_BMKDIR = join(BMK_DIR, "11-plot_hor_stv")


def get_hor_stv_plot_layout(wc) -> str:
    if config["humas_annot"]["mode"] == "sf":
        return workflow.source_path("../scripts/cenplot_sf_plot.yaml")
    elif config["humas_annot"]["mode"] == "srf-n-trf":
        return workflow.source_path("../scripts/cenplot_srf-n-trf_plot.yaml")
    else:
        return workflow.source_path("../scripts/cenplot_hor_stv_plot.yaml")


# TODO: Could be generalized to just be any pattern but lazy.
include: "11.1-plot_hor_stv_chr.smk"
include: "11.1-plot_hor_stv_sm.smk"


PLOT_HOR_STV_OUTPUTS = []
if "chromosome" in config["plot_hor_stv"].get("partition_by", []):
    PLOT_HOR_STV_OUTPUTS.extend(
        expand(
            rules.plot_hor_stv_by_chr.output,
            chr=CHROMOSOMES if CHROMOSOMES else ["all"],
            typ=["complete", "all"],
        )
    )
if "sample" in config["plot_hor_stv"].get("partition_by", []):
    PLOT_HOR_STV_OUTPUTS.extend(
        expand(
            rules.plot_hor_stv_by_sm.output,
            sm=SAMPLE_NAMES,
            typ=["complete", "all"],
        )
    )


rule plot_hor_stv_all:
    input:
        PLOT_HOR_STV_OUTPUTS,
    default_target: True
