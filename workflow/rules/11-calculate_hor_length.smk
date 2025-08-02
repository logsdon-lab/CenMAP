include: "common.smk"
include: "10-format_hor_stv.smk"


HOR_ARR_LEN_OUTDIR = join(OUTPUT_DIR, "11-calculate_hor_length")
HOR_ARR_LEN_LOGDIR = join(LOG_DIR, "11-calculate_hor_length")


rule calculate_as_hor_length:
    input:
        stv_row_bed=rules.filter_complete_correct_stv_row.output,
    output:
        stv_row_live_bed=join(
            HOR_ARR_LEN_OUTDIR,
            "bed",
            "{chr}_AS-HOR_stv_row.live.bed",
        ),
        arr_lens_bed=join(
            HOR_ARR_LEN_OUTDIR,
            "bed",
            "{chr}_AS-HOR_lengths.bed",
        ),
    conda:
        "../envs/py.yaml"
    log:
        join(HOR_ARR_LEN_LOGDIR, "calculate_{chr}_as_hor_length.log"),
    shell:
        """
        # Save live output.
        awk -v OFS="\\t" '$4 ~ "L"' {input.stv_row_bed} > {output.stv_row_live_bed}
        {{ censtats length -i {output.stv_row_live_bed} || true ;}} > {output.arr_lens_bed} 2> {log}
        """


rule aggregate_as_hor_length:
    input:
        expand(rules.calculate_as_hor_length.output.arr_lens_bed, chr=CHROMOSOMES),
    output:
        join(
            FINAL_OUTPUT_DIR,
            "bed",
            "all_AS-HOR_lengths.bed",
        ),
    shell:
        """
        cat {input} > {output}
        """


def args_ref_hor_lengths(wc, input) -> str:
    """
    Construct args for plotting ref hor length.
    """
    if not input.ref_hor_lengths:
        return ""
    args = []
    for ref in config["calculate_hor_length"]["ref_hor_lengths"]:
        lbl = ref["name"]
        path = ref["path"]
        color = ref["color"]
        # label, path, and color of violin plot dot.
        arg = f"{lbl}={path}={color}" if color else f"{lbl}={path}"
        args.append(arg)
    return f'-a {" ".join(args)}'


rule plot_as_hor_length:
    input:
        ref_hor_lengths=[
            ref["path"]
            for ref in config["calculate_hor_length"].get("ref_hor_lengths", [])
        ],
        lengths=rules.aggregate_as_hor_length.output,
    output:
        join(
            HOR_ARR_LEN_OUTDIR,
            "plots",
            "all_AS-HOR_lengths.png",
        ),
    params:
        script=workflow.source_path("../scripts/plot_hor_length.py"),
        # Random colors given for each ref.
        # Run separately to get desired colors.
        args_added_lengths=args_ref_hor_lengths,
    log:
        join(HOR_ARR_LEN_LOGDIR, "plot_all_as_hor_length.log"),
    conda:
        "../envs/py.yaml"
    shell:
        """
        if [ -s {input.lengths} ]; then
            python {params.script} \
            -i {input.lengths} \
            -o {output} {params.args_added_lengths} 2> {log}
        else
            touch {output}
        fi
        """


rule calculate_as_hor_length_all:
    input:
        rules.aggregate_as_hor_length.output,
        rules.plot_as_hor_length.output,
    default_target: True
