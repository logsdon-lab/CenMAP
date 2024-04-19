include: "common.smk"


rule calculate_as_hor_length:
    input:
        script="workflow/scripts/calculate_HOR_length.py",
        fmt_hmmer_output=os.path.join(
            config["plot_hor_stv"]["output_dir"], "bed", "{chr}_AS-HOR_stv_row.all.bed"
        ),
    output:
        os.path.join(
            config["calculate_hor_length"]["output_dir"],
            "lengths",
            "{chr}_AS-HOR_lengths.tsv",
        ),
    conda:
        "../env/py.yaml"
    log:
        "logs/calculate_{chr}_as_hor_length.log",
    params:
        bp_jump_thr=100_000,
        arr_len_thr=30_000,
    shell:
        """
        ( python {input.script} \
        --input {input.fmt_hmmer_output} \
        --bp_jump_thr {params.bp_jump_thr} \
        --arr_len_thr {params.arr_len_thr} || true ) > {output} 2> {log}
        """


rule aggregate_as_hor_length:
    input:
        expand(rules.calculate_as_hor_length.output, chr=CHROMOSOMES),
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
        script="workflow/scripts/plot_HOR_length.R",
        chm1_lengths=config["calculate_hor_length"]["chm1_hor_lengths"],
        chm13_lengths=config["calculate_hor_length"]["chm13_hor_lengths"],
        lengths=rules.aggregate_as_hor_length.output,
    output:
        os.path.join(
            config["calculate_hor_length"]["output_dir"],
            "plots",
            "all_AS-HOR_lengths.png",
        ),
    log:
        "logs/plot_all_as_hor_length.log",
    conda:
        "../env/r.yaml"
    shell:
        """
        Rscript {input.script} \
        --input {input.lengths} \
        --input_chm1 {input.chm1_lengths} \
        --input_chm13 {input.chm13_lengths} \
        --output {output} 2> {log}
        """


rule get_hor_length_only:
    input:
        rules.aggregate_as_hor_length.output,
        rules.plot_as_hor_length.output,
