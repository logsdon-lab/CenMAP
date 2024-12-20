
include: "common.smk"
include: "plot_hor_stv.smk"


if config.get("cdr_finder"):

    include: "cdr_finder.smk"


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
        chr_stv_row_bed=rules.aggregate_format_all_stv_row.output,
        chr_stv_row_ort_bed=rules.get_stv_row_ort_bed.output,
        all_sat_annot_bed=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "bed",
            "all_cens.annotation.bed",
        ),
        all_cdr_bed=lambda wc: (
            rules.merge_cdr_beds.output.reorient_cdr_output
            if config.get("cdr_finder")
            else []
        ),
        all_binned_methyl_bed=lambda wc: (
            rules.merge_cdr_beds.output.reorient_methyl_cdr_output
            if config.get("cdr_finder")
            else []
        ),
    output:
        sat_annot_bed=temp(
            os.path.join(OUTPUT_MODDOTPLOT_DIR, "{chr}_{fname}_sat_annot.bed")
        ),
        stv_row_bed=temp(
            os.path.join(
                OUTPUT_MODDOTPLOT_DIR,
                "{chr}_{fname}_stv_row.bed",
            )
        ),
        stv_row_ort_bed=temp(
            os.path.join(
                OUTPUT_MODDOTPLOT_DIR,
                "{chr}_{fname}_stv_row_ort.bed",
            )
        ),
        cdr_bed=temp(
            os.path.join(
                OUTPUT_MODDOTPLOT_DIR,
                "{chr}_{fname}_cdrs.bed",
            )
        ),
        binned_methyl_bed=temp(
            os.path.join(
                OUTPUT_MODDOTPLOT_DIR,
                "{chr}_{fname}_methyl.bed",
            )
        ),
    params:
        cdr_output=bool(config.get("cdr_finder", False)),
        # Use base of fname incase of off by ones in name coords
        fname_base=lambda wc: wc.fname.split(":")[0],
    shell:
        """
        ( grep '{params.fname_base}' {input.all_sat_annot_bed} || true ) > {output.sat_annot_bed}
        ( grep '{params.fname_base}' {input.chr_stv_row_bed} || true ) > {output.stv_row_bed}
        ( grep '{params.fname_base}' {input.chr_stv_row_ort_bed} || true) > {output.stv_row_ort_bed}
        if [ "{params.cdr_output}" != "False" ]; then
            ( grep '{params.fname_base}' {input.all_cdr_bed} || true ) > {output.cdr_bed}
            ( grep '{params.fname_base}' {input.all_binned_methyl_bed} || true ) > {output.binned_methyl_bed}
        else
            touch {output.cdr_bed}
            touch {output.binned_methyl_bed}
        fi
        """


rule plot_cen_moddotplot:
    input:
        script="workflow/scripts/plot_cens_moddotplot.R",
        seq_ident_bed=rules.run_moddotplot.output.bed,
        sat_annot_bed=rules.filter_annotations_moddotplot.output.sat_annot_bed,
        stv_row_bed=rules.filter_annotations_moddotplot.output.stv_row_bed,
        cdr_bed=rules.filter_annotations_moddotplot.output.cdr_bed,
        stv_row_ort_bed=rules.filter_annotations_moddotplot.output.stv_row_ort_bed,
        binned_methyl_bed=rules.filter_annotations_moddotplot.output.binned_methyl_bed,
    output:
        plots=expand(
            os.path.join(
                OUTPUT_MODDOTPLOT_DIR,
                "combined",
                "{{chr}}",
                "{{fname}}.tri.{ext}",
            ),
            ext=["png", "pdf"],
        ),
    params:
        mer_order=lambda wc: MONOMER_ORDER[wc.chr],
        output_dir=lambda wc, output: os.path.dirname(output.plots[0]),
    conda:
        "../envs/r.yaml"
    log:
        "logs/plot_cen_moddotplot/plot_cen_moddotplot_{chr}_{fname}.log",
    shell:
        """
        Rscript {input.script} \
        --bed {input.seq_ident_bed} \
        --hor {input.stv_row_bed} \
        --hor_ort {input.stv_row_ort_bed} \
        --sat {input.sat_annot_bed} \
        --cdr {input.cdr_bed} \
        --methyl {input.binned_methyl_bed} \
        --mer_order {params.mer_order} \
        --outdir {params.output_dir} 2>> {log}
        """


# https://stackoverflow.com/a/63040288
def moddotplot_outputs(wc):
    if config["moddotplot"].get("input_dir") is None:
        # Wait until done.
        _ = checkpoints.aggregate_format_all_stv_row.get(**wc).output

        fnames, chrs = extract_fnames_and_chr(
            os.path.join(config["humas_hmmer"]["input_dir"], "{fname}.fa"),
            filter_chr=wc.chr,
        )
    else:
        fnames, chrs = extract_fnames_and_chr(
            os.path.join(INPUT_FA_DIR, "{fname}.fa"),
            filter_chr=wc.chr,
        )

    return dict(
        moddotplot=expand(rules.run_moddotplot.output, zip, chr=chrs, fname=fnames),
        cen_moddoplot=expand(
            expand(
                rules.plot_cen_moddotplot.output,
                zip,
                fname=fnames,
                chr=chrs,
            ),
        ),
    )


rule moddotplot_all:
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


rule moddotplot_only:
    input:
        expand(rules.moddotplot_all.output, chr=CHROMOSOMES),
        rules._force_moddotplot_env_inclusion.output if IS_CONTAINERIZE_CMD else [],
    default_target: True
