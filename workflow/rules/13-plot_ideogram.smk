include: "common.smk"
include: "2-concat_asm.smk"
include: "9-get_complete_correct_cens.smk"


IDEOGRAM_OUTDIR = join(OUTPUT_DIR, "13-plot_ideogram")
IDEOGRAM_LOGDIR = join(LOG_DIR, "13-plot_ideogram")
IDEOGRAM_BMKDIR = join(BMK_DIR, "13-plot_ideogram")


def cmd_complete_correct_bed(wc, input) -> str:
    # Already formatted.
    if wc.typ == "complete":
        return input.bed

    # (name, st, end, ctg_name, ctg_len)
    cmd = [
        "awk",
        "-v",
        "OFS='\\t'",
        """{
            adj_name=$1;
            adj_st=$2;
            adj_end=$3;
            ctg_name=$4;
            ctg_len=$5;
            ort=($1 ~ "rc-") ? "-" : "+";
            ctg_st=$2; ctg_end=$3;
            if (ort == "-") {
                new_ctg_st=ctg_len-ctg_end + 1;
                new_ctg_end=ctg_len-ctg_st + 1;
                ctg_st=new_ctg_st;
                ctg_end=new_ctg_end;
            };
            print ctg_name, ctg_st, ctg_end, adj_name, 0, ort, adj_st, adj_end, "0,0,0"
        }""",
        input.bed,
    ]
    return f"<({' '.join(cmd)})"


rule plot_ideogram:
    input:
        bed=lambda wc: (
            rules.get_complete_correct_cens_bed.output.bed
            if wc.typ == "complete"
            else rules.make_complete_cens_bed.output.cen_bed
        ),
        fai=rules.concat_asm.output.idx,
    output:
        join(
            IDEOGRAM_OUTDIR,
            "{sm}_{typ}_ideogram.pdf",
        ),
    params:
        script=workflow.source_path("../scripts/plot_ideogram.py"),
        height=config["ideogram"]["height"],
        width=config["ideogram"]["width"],
        title=config["ideogram"]["title"],
        legend_prop=config["ideogram"]["legend_prop"],
        bed=cmd_complete_correct_bed,
        use_renamed_reoriented=(
            "--use_renamed_reoriented"
            if config["ideogram"]["use_renamed_reoriented"]
            else ""
        ),
    log:
        join(IDEOGRAM_LOGDIR, "plot_{sm}_{typ}_ideogram.log"),
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {params.script} \
        -i {params.bed} \
        -g {input.fai} \
        -o {output} \
        -t {params.title} \
        -ht {params.height} \
        -w {params.width} \
        --legend_prop {params.legend_prop} \
        {params.use_renamed_reoriented} 2> {log}
        """


rule plot_ideogram_all:
    input:
        expand(
            rules.plot_ideogram.output,
            sm=SAMPLE_NAMES,
            typ=["complete", "all"],
        ),
    default_target: True
