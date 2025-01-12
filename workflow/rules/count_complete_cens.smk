include: "common.smk"


rule count_complete_cens:
    input:
        script="workflow/scripts/count_complete_cens.py",
        hor_arr_len=os.path.join(
            config["calculate_hor_length"]["output_dir"],
            "lengths",
            "all_AS-HOR_lengths.tsv",
        ),
    output:
        cmp_cnts=os.path.join(
            config["count_complete_cens"]["output_dir"],
            "complete_cen_counts.tsv",
        ),
    log:
        "logs/count_cens/count_complete_cens.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        {{ python {input.script} -i {input.hor_arr_len} || true ;}} > {output} 2> {log}
        """


rule plot_complete_cen_counts:
    input:
        script="workflow/scripts/plot_complete_cen_counts.R",
        cmp_cnts=rules.count_complete_cens.output,
    output:
        os.path.join(
            config["count_complete_cens"]["output_dir"],
            "complete_cen_counts.png",
        ),
    params:
        label=config["count_complete_cens"]["plot_lbl"],
        color=config["count_complete_cens"]["plot_color"],
    log:
        "logs/count_cens/plot_complete_cen_counts.log",
    conda:
        "../envs/r.yaml"
    shell:
        """
        if [ -s {input.cmp_cnts} ]; then
            Rscript {input.script} --input {input.cmp_cnts} --output {output} 2> {log}
        else
            touch {output}
        fi
        """


rule count_complete_cens_only:
    input:
        rules.count_complete_cens.output,
        rules.plot_complete_cen_counts.output,
