
include: "common.smk"
include: "utils.smk"
## TODO: When module rule inheritance is fixed, revert.
# include: "8-cdr_finder.smk"
include: "8-format_repeatmasker_sat_annot.smk"
include: "10-format_hor_stv.smk"


MODDOTPLOT_OUTDIR = join(OUTPUT_DIR, "11-moddotplot")
MODDOTPLOT_LOGDIR = join(LOG_DIR, "11-moddotplot")
MODDOTPLOT_BMKDIR = join(BMK_DIR, "11-moddotplot")
MODDOTPLOT_OUTDIR_OG = join(MODDOTPLOT_OUTDIR, "original")
MODDOTPLOT_OUTDIR_DATA = join(MODDOTPLOT_OUTDIR, "combined", "bed")
MODDOTPLOT_OUTDIR_PLOT = join(MODDOTPLOT_OUTDIR, "combined", "plots")


# Use ModDotPlot's ANI algorithm to assess self-identity.
# I would use ModDotPlot but no conda support, large dependencies, and pysam issues forces me to use my own implementation.
rule run_moddotplot:
    input:
        fasta=join(HUMAS_CENS_SPLIT_DIR, "{fname}.fa"),
    output:
        bedpe=join(MODDOTPLOT_OUTDIR_OG, "{chr}", "{fname}", "{fname}.bed"),
    conda:
        "../envs/py.yaml"
    params:
        window=config["moddotplot"]["window"],
        ident_thr=config["moddotplot"]["ident_thr"],
        outdir=lambda wc, output: os.path.dirname(output.bedpe),
    resources:
        mem=config["moddotplot"]["mem"],
    log:
        join(MODDOTPLOT_LOGDIR, "moddotplot_{chr}_{fname}.log"),
    benchmark:
        join(MODDOTPLOT_BMKDIR, "moddotplot_{chr}_{fname}.tsv")
    shell:
        """
        censtats self-ident -i {input.fasta} -x 2D -w {params.window} -o {params.outdir} -t {params.ident_thr} &> {log}
        """


rule filter_annotations_moddotplot:
    input:
        chr_stv_row_bed=rules.aggregate_format_all_stv_row.output,
        all_sat_annot_bed=rules.create_rm_satellite_annotations.output,
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
        all_binned_methyl_bed=(
            join(
                OUTPUT_DIR,
                "8-cdr_finder",
                "bed",
                "all_binned_freq.bed",
            )
            if config.get("cdr_finder")
            else []
        ),
        ident_bedpe=rules.run_moddotplot.output.bedpe,
    output:
        sat_annot_bed=join(MODDOTPLOT_OUTDIR_DATA, "{chr}", "{fname}", "sat_annot.bed"),
        ident_bedpe=join(MODDOTPLOT_OUTDIR_DATA, "{chr}", "{fname}", "ident.bedpe"),
        stv_row_bed=join(
            MODDOTPLOT_OUTDIR_DATA,
            "{chr}",
            "{fname}",
            "stv.bed",
        ),
        cdr_bed=join(
            MODDOTPLOT_OUTDIR_DATA,
            "{chr}",
            "{fname}",
            "cdrs.bed",
        ),
        binned_methyl_bed=join(
            MODDOTPLOT_OUTDIR_DATA,
            "{chr}",
            "{fname}",
            "methyl.bed",
        ),
    params:
        cdr_output=bool(config.get("cdr_finder", False)),
    shell:
        """
        ( grep '{wildcards.fname}' {input.all_sat_annot_bed} || true ) > {output.sat_annot_bed}
        ( grep '{wildcards.fname}' {input.chr_stv_row_bed} || true ) > {output.stv_row_bed}
        # Remove header.
        ln -s $(realpath {input.ident_bedpe}) {output.ident_bedpe}
        if [ "{params.cdr_output}" != "False" ]; then
            ( grep '{wildcards.fname}' {input.all_cdr_bed} || true ) > {output.cdr_bed}
            ( grep '{wildcards.fname}' {input.all_binned_methyl_bed} || true ) > {output.binned_methyl_bed}
        else
            touch {output.cdr_bed}
            touch {output.binned_methyl_bed}
        fi
        """


def get_moddotplot_plot_layout(wc) -> str:
    if config["humas_annot"]["mode"] == "sf":
        return workflow.source_path("../scripts/cenplot_moddotplot_sf.yaml")
    elif config["humas_annot"]["mode"] == "srf-n-trf":
        return workflow.source_path("../scripts/cenplot_moddotplot_srf-n-trf.yaml")
    else:
        return workflow.source_path("../scripts/cenplot_moddotplot.yaml")


rule modify_moddotplot_cenplot_tracks:
    input:
        **(
            {"cdr": rules.filter_annotations_moddotplot.output.cdr_bed}
            if config.get("cdr_finder")
            else {}
        ),
        rm=rules.filter_annotations_moddotplot.output.sat_annot_bed,
        ident=rules.filter_annotations_moddotplot.output.ident_bedpe,
        stv=rules.filter_annotations_moddotplot.output.stv_row_bed,
        methyl=rules.filter_annotations_moddotplot.output.binned_methyl_bed,
    output:
        plot_layout=join(
            MODDOTPLOT_OUTDIR_PLOT,
            "{chr}",
            "{fname}.yaml",
        ),
    conda:
        "../envs/py.yaml"
    log:
        join(MODDOTPLOT_LOGDIR, "modify_moddotplot_cenplot_tracks_{chr}_{fname}.log"),
    params:
        script=workflow.source_path("../scripts/format_cenplot_yaml.py"),
        plot_layout=get_moddotplot_plot_layout,
        json_file_str=lambda wc, input: json.dumps(dict(input)),
        options=lambda wc: f"--options '{json.dumps({"hor":{"color_map_file": os.path.abspath(config["plot_hor_stv"]["stv_annot_colors"])}})}'",
    shell:
        """
        python {params.script} \
        -i '{params.json_file_str}' \
        -t {params.plot_layout} \
        {params.options} > {output} 2> {log}
        """


rule plot_cen_moddotplot:
    input:
        plot_layout=rules.modify_moddotplot_cenplot_tracks.output,
        bedfiles=rules.filter_annotations_moddotplot.output,
    output:
        plots=multiext(
            join(
                MODDOTPLOT_OUTDIR_PLOT,
                "{chr}",
                "{fname}",
            ),
            ".pdf",
            ".png",
        ),
    params:
        output_dir=lambda wc, output: os.path.dirname(output.plots[0]),
    conda:
        "../envs/py.yaml"
    log:
        join(MODDOTPLOT_LOGDIR, "plot_moddotplot_{chr}_{fname}.log"),
    shell:
        """
        cenplot draw -t {input.plot_layout} -d {params.output_dir} -c {wildcards.fname} -p {threads} 2> {log}
        """


# https://stackoverflow.com/a/63040288
def moddotplot_outputs(wc):
    # Wait until done.
    try:
        _ = checkpoints.aggregate_format_all_stv_row.get(**wc).output
    except AttributeError:
        pass

    fastas = glob.glob(join(HUMAS_CENS_SPLIT_DIR, "*.fa"))
    fnames = get_valid_fnames(fastas, filter_chrom=wc.chr if wc.chr != "all" else None)

    return dict(
        moddotplot=expand(rules.run_moddotplot.output, chr=wc.chr, fname=fnames),
        cen_moddoplot=expand(
            expand(
                rules.plot_cen_moddotplot.output,
                fname=fnames,
                chr=wc.chr,
            ),
        ),
    )


rule moddotplot_chr:
    input:
        unpack(moddotplot_outputs),
    output:
        touch(join(MODDOTPLOT_OUTDIR, "moddotplot_{chr}.done")),


# Force moddotplot to be included with --containerize.
rule _force_moddotplot_env_inclusion:
    output:
        plots=touch("conda_moddotplot.done"),
    conda:
        "../envs/py.yaml"
    shell:
        "echo ''"


rule moddotplot_all:
    input:
        expand(
            rules.moddotplot_chr.output,
            chr=CHROMOSOMES if CHROMOSOMES else ["all"],
        ),
        rules._force_moddotplot_env_inclusion.output if IS_CONTAINERIZE_CMD else [],
    default_target: True
