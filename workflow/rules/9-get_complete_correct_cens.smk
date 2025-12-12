include: "common.smk"
include: "utils.smk"
include: "7-finalize_cens.smk"


COMPLETE_CORRECT_CENS_LOGDIR = join(LOG_DIR, "9-get_complete_correct_cens")


if config.get("nucflag"):

    include: "8-nucflag.smk"


if IS_HUMAN_ANNOT:

    include: "8-humas_annot.smk"


def cmd_infile(wc, input):
    cmd = []
    if input.nucflag_bed:
        cmd.extend(
            [
                "bedtools",
                "intersect",
                "-wa",
                "-u",
                "-a",
                f"<(sort -k1,1 -k2,2n {input.interm_bed})",
                "-b",
                f"<(sort -k1,1 -k2,2n {input.nucflag_bed} | awk '$4 == \"good\"')",
            ]
        )
    else:
        cmd.extend(["cat", input.interm_bed])

    if (
        config.get("get_complete_correct_cens", {}).get("filter_by_hor", True)
        and "cmd_intersect_live" in globals()
    ):
        pipe_intersect_live = cmd_intersect_live(wc)
        cmd.append(pipe_intersect_live)
    return " ".join(cmd)


rule get_complete_correct_cens_bed:
    input:
        # (name, st, end, ctg_name, ctg_len)
        interm_bed=rules.make_complete_cens_bed.output.cen_bed,
        # (name, st, end, status)
        nucflag_bed=(
            rules.check_asm_nucflag.output.asm_status if config.get("nucflag") else []
        ),
        stv_chkpt=rules.sm_stv.output if IS_HUMAN_ANNOT else [],
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
        infile_stream=cmd_infile,
        cmd_into_bed9=cmd_complete_bed_as_bed9(),
    log:
        join(
            COMPLETE_CORRECT_CENS_LOGDIR,
            "get_complete_correct_bed_{sm}.log",
        ),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ {params.infile_stream} | \
        {params.cmd_into_bed9} ;}} > {output} 2> {log}
        """


rule get_complete_correct_cens_fa:
    input:
        # assembly reoriented fasta
        fa=rules.create_final_asm.output.fa,
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
        bed=lambda wc, input: f"""<(awk -v OFS="\\t" '{{new_len=$7-1; print $4, (new_len < 0) ? 0 : new_len, $8}}' {input.bed})""",
        added_cmds=f"| {cmd_filter_fa_chrom()}" if CHROMOSOMES else "",
    shell:
        shell_extract_and_index_fa


rule get_complete_correct_cens_all:
    input:
        expand(rules.get_complete_correct_cens_bed.output, sm=SAMPLE_NAMES),
        expand(rules.get_complete_correct_cens_fa.output, sm=SAMPLE_NAMES),
    default_target: True
