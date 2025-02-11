include: "common.smk"


rule calculate_as_hor_length:
    input:
        stv_row_bed=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "bed",
            "{chr}_AS-HOR_stv_row.complete.bed",
        ),
    output:
        stv_row_live_bed=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "bed",
            "{chr}_AS-HOR_stv_row.live.bed",
        ),
        arr_lens_bed=os.path.join(
            config["calculate_hor_length"]["output_dir"],
            "lengths",
            "{chr}_AS-HOR_lengths.tsv",
        ),
    conda:
        "../envs/py.yaml"
    log:
        "logs/calculate_hor_length/calculate_{chr}_as_hor_length.log",
    shell:
        """
        awk -v OFS="\\t" '$4 ~ "L"' {input.stv_row_bed} > {output.stv_row_live_bed}
        {{ censtats length -i {output.stv_row_live_bed} || true ;}} > {output.arr_lens_bed} 2> {log}
        """


rule aggregate_as_hor_length:
    input:
        expand(rules.calculate_as_hor_length.output.arr_lens_bed, chr=CHROMOSOMES),
    output:
        os.path.join(
            config["calculate_hor_length"]["output_dir"],
            "lengths",
            "all_AS-HOR_lengths.tsv",
        ),
    shell:
        """
        cat {input} > {output}
        """


rule plot_as_hor_length:
    input:
        script="workflow/scripts/plot_hor_length.py",
        comparison_hor_lengths=config["calculate_hor_length"].get(
            "comparison_hor_lengths", []
        ),
        lengths=rules.aggregate_as_hor_length.output,
    output:
        os.path.join(
            config["calculate_hor_length"]["output_dir"],
            "plots",
            "all_AS-HOR_lengths.png",
        ),
    params:
        args_added_lengths=lambda wc, input: "-a "
        + " ".join(input.comparison_hor_lengths)
        if input.comparison_hor_lengths
        else "",
    log:
        "logs/calculate_hor_length/plot_all_as_hor_length.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        if [ -s {input.lengths} ]; then
            python {input.script} \
            -i {input.lengths} \
            -o {output} {params.comparison_hor_lengths} 2> {log}
        else
            touch {output}
        fi
        """


rule calculate_as_hor_length_all:
    input:
        rules.aggregate_as_hor_length.output,
        rules.plot_as_hor_length.output,
