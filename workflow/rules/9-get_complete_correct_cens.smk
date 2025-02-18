include: "common.smk"
include: "utils.smk"
include: "7-fix_cens_w_repeatmasker.smk"


COMPLETE_CORRECT_CENS_LOGDIR = join(LOG_DIR, "9-get_complete_correct_cens")


def nucflag_bed(wc):
    if not config.get("nucflag"):
        return []

    # Include only if nucflag added.
    include: "8-nucflag.smk"

    return rules.check_asm_nucflag.output.asm_status


rule get_complete_correct_cens_bed:
    input:
        # (name, st, end, is_partial, ctg_name, ctg_len)
        interm_bed=rules.make_complete_cens_bed.output,
        # (name, st, end, status)
        nucflag_bed=nucflag_bed,
    output:
        # BED9
        # (name, st, end, adj_name, score, ort, adj_st, adj_end, rgb)
        bed=join(
            FINAL_OUTPUT_DIR,
            "bed",
            "{sm}_complete_correct_cens.bed",
        ),
    params:
        # Allow not running nucflag.
        infile_stream=lambda wc, input: (
            f"join <(sort -k1 {input.interm_bed}) <(sort -k1 {input.nucflag_bed})"
            if input.nucflag_bed
            else f"cat {input.interm_bed}"
        ),
    log:
        join(COMPLETE_CORRECT_CENS_LOGDIR, "get_complete_correct_bed_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ {params.infile_stream} | \
        awk -v OFS="\\t" '{{
            adj_name=$1;
            adj_st=$2; adj_end=$3;
            ctg_name=$5;
            ctg_len=$6;
            ort=($1 ~ "rc-") ? "-" : "+";
            ctg_st=$2; ctg_end=$3;
            if (ort == "-") {{
                new_ctg_st=ctg_len-ctg_end + 1;
                new_ctg_end=ctg_len-ctg_st + 1;
                ctg_st=new_ctg_st;
                ctg_end=new_ctg_end;
            }}
            status=(($9 == "good" || $9 == "") && ($4 == "false")) ? "good" : "misassembled";
            if (status != "good") {{
                next;
            }};
            print ctg_name, ctg_st, ctg_end, adj_name, 0, ort, adj_st, adj_end, "0,0,0"
        }}' ;}} > {output} 2> {log}
        """


rule get_complete_correct_cens_fa:
    input:
        # assembly reoriented fasta
        fa=rules.rename_reort_asm.output.fa,
        # (name, st, end)
        bed=rules.get_complete_correct_cens_bed.output,
    output:
        seq=join(
            FINAL_OUTPUT_DIR,
            "seq",
            "{sm}_cens.fa",
        ),
        idx=join(
            FINAL_OUTPUT_DIR,
            "seq",
            "{sm}_cens.fa.fai",
        ),
    conda:
        "../envs/tools.yaml"
    log:
        join(COMPLETE_CORRECT_CENS_LOGDIR, "get_complete_correct_cens_fa_{sm}.log"),
    params:
        added_cmds="",
        bed=lambda wc, input: f"""<(awk -v OFS="\\t" '{{print $4, $7-1, $8}}' {input.bed})""",
    shell:
        shell_extract_and_index_fa


rule get_complete_correct_cens_all:
    input:
        expand(rules.get_complete_correct_cens_bed.output, sm=SAMPLE_NAMES),
        expand(rules.get_complete_correct_cens_fa.output, sm=SAMPLE_NAMES),
    default_target: True
