
include: "common.smk"
include: "utils.smk"


if config["moddotplot"].get("input_dir"):
    INPUT_FA_DIR = config["moddotplot"]["input_dir"]
else:
    INPUT_FA_DIR = HUMAS_CENS_SPLIT_DIR

OUTPUT_MODDOTPLOT_DIR = config["moddotplot"].get("output_dir", "results/moddotplot")
OUTPUT_MODDOTPLOT_DATA_DIR = os.path.join(OUTPUT_MODDOTPLOT_DIR, "combined", "bed")
OUTPUT_MODDOTPLOT_PLOTS_DIR = os.path.join(OUTPUT_MODDOTPLOT_DIR, "combined", "plots")


rule run_moddotplot:
    input:
        fasta=os.path.join(INPUT_FA_DIR, "{fname}.fa"),
    output:
        plots=expand(
            os.path.join(
                OUTPUT_MODDOTPLOT_DIR,
                "original",
                "{{chr}}",
                "{{fname}}",
                "{{fname}}_{otype}.{ext}",
            ),
            otype=["FULL", "HIST", "TRI"],
            ext=["png", "pdf"],
        ),
        bed=os.path.join(
            OUTPUT_MODDOTPLOT_DIR, "original", "{chr}", "{fname}", "{fname}.bed"
        ),
    conda:
        "../envs/py.yaml"
    params:
        window=config["moddotplot"]["window"],
        ident_thr=70.0,
        outdir=lambda wc, output: os.path.dirname(output.bed),
    resources:
        mem=config["moddotplot"]["mem"],
    log:
        "logs/moddotplot/moddotplot_{chr}_{fname}.log",
    benchmark:
        "benchmarks/moddotplot/moddotplot_{chr}_{fname}.tsv"
    shell:
        """
        moddotplot static -f {input.fasta} -w {params.window} -o {params.outdir} -id {params.ident_thr} &> {log}
        """


rule filter_annotations_moddotplot:
    input:
        chr_stv_row_bed=os.path.join(
            config["plot_hor_stv"]["output_dir"], "bed", "{chr}", "stv_all.bed"
        ),
        all_sat_annot_bed=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "bed",
            "all_cens.annotation.bed",
        ),
        all_cdr_bed=lambda wc: (
            os.path.join(
                config["cdr_finder"]["output_dir"],
                "bed",
                "all_cdrs.bed",
            )
            if config.get("cdr_finder")
            else []
        ),
        all_binned_methyl_bed=lambda wc: (
            os.path.join(
                config["cdr_finder"]["output_dir"],
                "bed",
                "all_binned_freq.bed",
            )
            if config.get("cdr_finder")
            else []
        ),
        ident_bed=rules.run_moddotplot.output.bed,
    output:
        sat_annot_bed=os.path.join(
            OUTPUT_MODDOTPLOT_DATA_DIR, "{chr}", "{fname}", "sat_annot.bed"
        ),
        ident_bed=os.path.join(
            OUTPUT_MODDOTPLOT_DATA_DIR, "{chr}", "{fname}", "ident.bed"
        ),
        stv_row_bed=os.path.join(
            OUTPUT_MODDOTPLOT_DATA_DIR,
            "{chr}",
            "{fname}",
            "stv.bed",
        ),
        cdr_bed=os.path.join(
            OUTPUT_MODDOTPLOT_DATA_DIR,
            "{chr}",
            "{fname}",
            "cdrs.bed",
        ),
        binned_methyl_bed=os.path.join(
            OUTPUT_MODDOTPLOT_DATA_DIR,
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


use rule modify_cenplot_tracks as modify_moddotplot_cenplot_tracks with:
    input:
        plot_layout="workflow/scripts/cenplot_moddotplot.toml",
        infiles=[rules.filter_annotations_moddotplot.output.cdr_bed],
    output:
        plot_layout=os.path.join(
            OUTPUT_MODDOTPLOT_PLOTS_DIR,
            "{chr}",
            "{fname}.yaml",
        ),


rule plot_cen_moddotplot:
    input:
        plot_layout=rules.modify_moddotplot_cenplot_tracks.output,
        bedfiles=rules.filter_annotations_moddotplot.output,
    output:
        plots=multiext(
            os.path.join(
                OUTPUT_MODDOTPLOT_PLOTS_DIR,
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
        "logs/moddotplot/plot_moddotplot_{chr}_{fname}.log",
    shell:
        """
        cenplot draw -t {input.plot_layout} -d {params.output_dir} -c <(echo {wildcards.fname}) -p {threads} 2> {log}
        """


# https://stackoverflow.com/a/63040288
def moddotplot_outputs(wc):
    if config["moddotplot"].get("input_dir") is None:
        # Wait until done.
        try:
            _ = checkpoints.aggregate_format_all_stv_row.get(**wc).output
        except AttributeError:
            pass

        wcs = glob_wildcards(
            os.path.join(HUMAS_CENS_SPLIT_DIR, "{sm}_" + wc.chr + "_{ctg}.fa"),
        )
    else:
        wcs = glob_wildcards(
            os.path.join(INPUT_FA_DIR, "{sm}_" + wc.chr + "_{ctg}.fa"),
        )
    fnames = [f"{sm}_{wc.chr}_{ctg}" for sm, ctg in zip(wcs.sm, wcs.ctg)]
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
        touch(os.path.join(OUTPUT_MODDOTPLOT_DIR, "moddotplot_{chr}.done")),


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
