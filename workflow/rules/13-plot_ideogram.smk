include: "common.smk"
include: "2-concat_asm.smk"
include: "9-get_complete_correct_cens.smk"


IDEOGRAM_OUTDIR = join(OUTPUT_DIR, "13-plot_ideogram")
IDEOGRAM_LOGDIR = join(LOG_DIR, "13-plot_ideogram")
IDEOGRAM_BMKDIR = join(BMK_DIR, "13-plot_ideogram")


rule plot_ideogram:
    input:
        bed=rules.get_complete_correct_cens_bed.output.bed,
        fai=rules.concat_asm.output.idx,
    output:
        join(
            IDEOGRAM_OUTDIR,
            "{sm}_ideogram.pdf",
        ),
    params:
        script=workflow.source_path("../scripts/plot_ideogram.py"),
        height=config["ideogram"]["height"],
        width=config["ideogram"]["width"],
        title=config["ideogram"]["title"],
        legend_prop=config["ideogram"]["legend_prop"],
    log:
        join(IDEOGRAM_LOGDIR, "plot_{sm}_ideogram.log"),
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {params.script} \
        -i {input.bed} \
        -g {input.fai} \
        -o {output} \
        -t {params.title} \
        -ht {params.height} \
        -w {params.width} \
        --legend_prop {params.legend_prop} 2> {log}
        """


rule plot_ideogram_all:
    input:
        expand(rules.plot_ideogram.output, sm=SAMPLE_NAMES),
    default_target: True
