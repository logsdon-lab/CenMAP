
include: "common.smk"


if config["moddotplot"].get("input_dir"):
    INPUT_FA_DIR = config["moddotplot"]["input_dir"]
else:
    INPUT_FA_DIR = config["humas_hmmer"]["input_dir"]

OUTPUT_MODDOTPLOT_DIR = config["moddotplot"].get("output_dir", "results/moddotplot")


rule run_moddotplot:
    input:
        fasta=os.path.join(INPUT_FA_DIR, "{fname}.fa"),
    output:
        plots=expand(
            os.path.join(
                OUTPUT_MODDOTPLOT_DIR, "original_{{fname}}", "{{fname}}_{otype}.{ext}"
            ),
            otype=["FULL", "HIST", "TRI"],
            ext=["png", "pdf"],
        ),
        bed=os.path.join(OUTPUT_MODDOTPLOT_DIR, "original_{fname}", "{fname}.bed"),
    conda:
        "../env/moddotplot.yaml"
    params:
        window=config["moddotplot"]["window"],
        ident_thr=70.0,
        outdir=lambda wc: os.path.join(OUTPUT_MODDOTPLOT_DIR, f"original_{wc.fname}"),
    resources:
        mem=config["moddotplot"]["mem"],
    log:
        "logs/moddotplot/moddotplot_{fname}.log",
    benchmark:
        "benchmarks/moddotplot/moddotplot_{fname}.tsv"
    shell:
        """
        moddotplot static -f {input.fasta} -w {params.window} -o {params.outdir} -id {params.ident_thr} &> {log}
        """


rule plot_cen_moddotplot:
    input:
        script="workflow/scripts/repeats_moddotplot.R",
        seq_ident_bed=rules.run_moddotplot.output.bed,
        # TODO: Remove grepping. Pass file directly. Will require stv wf refactor.
        chr_stv_row_bed=os.path.join(
            config["plot_hor_stv"]["output_dir"], "bed", "{chr}_AS-HOR_stv_row.all.bed"
        ),
        all_sat_annot_bed=os.path.join(
            config["repeatmasker_sat_annot"]["output_dir"],
            "bed",
            "all_cens.annotation.bed",
        ),
    output:
        sat_annot_bed=temp(
            os.path.join(
                OUTPUT_MODDOTPLOT_DIR, "{chr}_{mer_order}_{fname}_sat_annot.bed"
            )
        ),
        stv_row_bed=temp(
            os.path.join(
                OUTPUT_MODDOTPLOT_DIR,
                "{chr}_{mer_order}_{fname}_stv_row.bed",
            )
        ),
        plots=expand(
            os.path.join(
                OUTPUT_MODDOTPLOT_DIR,
                "{{chr}}_{{fname}}_full",
                "{{fname}}_{{mer_order}}.tri.{ext}",
            ),
            ext=["png", "pdf"],
        ),
    params:
        output_dir=lambda wc: os.path.join(
            OUTPUT_MODDOTPLOT_DIR, f"{wc.chr}_{wc.fname}_full"
        ),
    conda:
        "../env/r.yaml"
    log:
        "logs/plot_cen_moddotplot/plot_cen_moddotplot_{chr}_{fname}_{mer_order}.log",
    shell:
        """
        grep '{wildcards.fname}' {input.all_sat_annot_bed} > {output.sat_annot_bed}
        grep '{wildcards.fname}' {input.chr_stv_row_bed} > {output.stv_row_bed}
        Rscript {input.script} \
        --bed {input.seq_ident_bed} \
        --hor {output.stv_row_bed} \
        --sat {output.sat_annot_bed} \
        --mer_order {wildcards.mer_order} \
        --outdir {params.output_dir} 2>> {log}
        """


# https://stackoverflow.com/a/63040288
def moddotplot_outputs_no_input_dir(wc):
    # Wait until done.
    try:
        _ = checkpoints.split_cens_for_humas_hmmer.get(**wc).output
    except AttributeError:
        pass

    fnames, chrs = extract_fa_fnames_and_chr(config["humas_hmmer"]["input_dir"])

    wildcard_constraints:
        fname="|".join(fnames),

    return dict(
        modotplot=expand(rules.run_moddotplot.output, fname=fnames),
        cen_moddoplot=expand(
            expand(
                rules.plot_cen_moddotplot.output,
                zip,
                fname=fnames,
                chr=chrs,
                mer_order=["{mer_order}"] * len(fnames),
            ),
            mer_order=config["plot_hor_stv"]["mer_order"],
        ),
    )


# Conditionally change based on provided input dir.
if config["moddotplot"].get("input_dir") is None:

    rule moddotplot_all:
        input:
            unpack(moddotplot_outputs_no_input_dir),
        output:
            touch(os.path.join(OUTPUT_MODDOTPLOT_DIR, "moddotplot_{chr}.done")),

else:
    # Extract filenames and chromosome names from input directory.
    fnames, chrs = extract_fa_fnames_and_chr(INPUT_FA_DIR)

    wildcard_constraints:
        fname="|".join(fnames),

    rule moddotplot_all:
        input:
            expand(
                rules.run_moddotplot.output,
                fname=fnames,
            ),
            expand(
                expand(
                    rules.plot_cen_moddotplot.output,
                    zip,
                    fname=fnames,
                    chr=chrs,
                    mer_order=["{mer_order}"] * len(fnames),
                ),
                mer_order=config["plot_hor_stv"]["mer_order"],
            ),


rule moddotplot_only:
    input:
        (
            expand(rules.moddotplot_all.output, chr=CHROMOSOMES)
            if config["moddotplot"].get("input_dir") is None
            else rules.moddotplot_all.input
        ),
    default_target: True
