
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


rule run_moddotplot:
    input:
        fasta=join(HUMAS_CENS_SPLIT_DIR, "{fname}.fa"),
    output:
        plots=expand(
            join(
                MODDOTPLOT_OUTDIR_OG,
                "{{chr}}",
                "{{fname}}",
                "{{fname}}_{otype}.{ext}",
            ),
            otype=["FULL", "HIST", "TRI"],
            ext=["png", "pdf"],
        ),
        bed=join(MODDOTPLOT_OUTDIR_OG, "{chr}", "{fname}", "{fname}.bed"),
    conda:
        "../envs/py.yaml"
    params:
        window=config["moddotplot"]["window"],
        ident_thr=config["moddotplot"]["ident_thr"],
        outdir=lambda wc, output: os.path.dirname(output.bed),
    resources:
        mem=config["moddotplot"]["mem"],
    log:
        join(MODDOTPLOT_LOGDIR, "moddotplot_{chr}_{fname}.log"),
    benchmark:
        join(MODDOTPLOT_BMKDIR, "moddotplot_{chr}_{fname}.tsv")
    shell:
        """
        moddotplot static -f {input.fasta} -w {params.window} -o {params.outdir} -id {params.ident_thr} &> {log}
        """


rule filter_annotations_moddotplot:
    input:
        chr_stv_row_bed=rules.aggregate_format_all_stv_row.output,
        all_sat_annot_bed=rules.aggregate_rm_satellite_annotations.output,
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
        ident_bed=rules.run_moddotplot.output.bed,
    output:
        sat_annot_bed=join(MODDOTPLOT_OUTDIR_DATA, "{chr}", "{fname}", "sat_annot.bed"),
        ident_bed=join(MODDOTPLOT_OUTDIR_DATA, "{chr}", "{fname}", "ident.bed"),
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
        # Use base of fname incase of off by ones in name coords
        fname_base=lambda wc: wc.fname.split(":")[0],
    shell:
        """
        ( grep '{params.fname_base}' {input.all_sat_annot_bed} || true ) > {output.sat_annot_bed}
        ( grep '{params.fname_base}' {input.chr_stv_row_bed} || true ) > {output.stv_row_bed}
        # Remove header.
        tail -n+2 {input.ident_bed} > {output.ident_bed}
        if [ "{params.cdr_output}" != "False" ]; then
            ( grep '{params.fname_base}' {input.all_cdr_bed} || true ) > {output.cdr_bed}
            ( grep '{params.fname_base}' {input.all_binned_methyl_bed} || true ) > {output.binned_methyl_bed}
        else
            touch {output.cdr_bed}
            touch {output.binned_methyl_bed}
        fi
        """


rule modify_moddotplot_cenplot_tracks:
    input:
        plot_layout=workflow.source_path("../scripts/cenplot_moddotplot.toml"),
        infiles=[rules.filter_annotations_moddotplot.output.cdr_bed],
    output:
        plot_layout=join(
            MODDOTPLOT_OUTDIR_PLOT,
            "{chr}",
            "{fname}.yaml",
        ),
    run:
        format_toml_path(
            input_plot_layout=input.plot_layout,
            output_plot_layout=output.plot_layout,
            indir=os.path.abspath(os.path.dirname(str(input.infiles[0]))),
            **dict(params.items()),
        )


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
    threads: 1
    params:
        output_dir=lambda wc, output: os.path.dirname(output.plots[0]),
    conda:
        "../envs/py.yaml"
    log:
        join(MODDOTPLOT_LOGDIR, "plot_moddotplot_{chr}_{fname}.log"),
    shell:
        """
        cenplot draw -t {input.plot_layout} -d {params.output_dir} -c <(echo {wildcards.fname}) -p {threads} 2> {log}
        """


# https://stackoverflow.com/a/63040288
def moddotplot_outputs(wc):
    # Wait until done.
    try:
        _ = checkpoints.aggregate_format_all_stv_row.get(**wc).output
    except AttributeError:
        pass

    wcs = glob_wildcards(
        join(HUMAS_CENS_SPLIT_DIR, "{sm}_{chr}_{ctg}.fa"),
    )
    fnames = [
        f"{sm}_{chrom}_{ctg}"
        for sm, chrom, ctg in zip(wcs.sm, wcs.chr, wcs.ctg)
        if chrom == wc.chr or chrom == f"rc-{wc.chr}"
    ]
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
        expand(rules.moddotplot_chr.output, chr=CHROMOSOMES),
        rules._force_moddotplot_env_inclusion.output if IS_CONTAINERIZE_CMD else [],
    default_target: True
