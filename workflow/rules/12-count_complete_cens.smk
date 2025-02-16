include: "common.smk"
include: "11-calculate_hor_length.smk"


CNT_CENS_OUTDIR = join(OUTPUT_DIR, "12-count_complete_cens")
CNT_CENS_LOGDIR = join(LOG_DIR, "12-count_complete_cens")
CNT_CENS_BMKDIR = join(BMK_DIR, "12-count_complete_cens")


rule count_complete_cens:
    input:
        script=workflow.source_path("../scripts/count_complete_cens.py"),
        hor_arr_len=rules.aggregate_as_hor_length.output,
    output:
        cmp_cnts=join(
            CNT_CENS_OUTDIR,
            "complete_cen_counts.tsv",
        ),
    log:
        join(CNT_CENS_LOGDIR, "count_complete_cens.log"),
    conda:
        "../envs/py.yaml"
    shell:
        """
        {{ python {input.script} -i {input.hor_arr_len} || true ;}} > {output} 2> {log}
        """


rule plot_complete_cen_counts:
    input:
        script=workflow.source_path("../scripts/plot_complete_cen_counts.py"),
        cmp_cnts=rules.count_complete_cens.output,
    output:
        join(
            CNT_CENS_OUTDIR,
            "complete_cen_counts.png",
        ),
    log:
        join(CNT_CENS_LOGDIR, "plot_complete_cen_counts.log"),
    conda:
        "../envs/py.yaml"
    shell:
        """
        if [ -s {input.cmp_cnts} ]; then
            python {input.script} -i {input.cmp_cnts} -o {output} 2> {log}
        else
            touch {output}
        fi
        """


rule count_complete_cens_all:
    input:
        rules.count_complete_cens.output,
        rules.plot_complete_cen_counts.output,
    default_target: True
