include: "nucflag.smk"


rule get_complete_correct_cens_bed:
    input:
        # (name_w_coords, st, end, is_partial, name)
        interm_bed=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            "interm",
            "{sm}_complete_cens_w_coords.bed",
        ),
        # (name_w_coords, st, end, status)
        nucflag_bed=(
            rules.check_asm_nucflag.output.asm_status if config.get("nucflag") else []
        ),
    output:
        # (name_no_coords, st, end)
        bed=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            "final",
            "{sm}_complete_correct_cens.bed",
        ),
    params:
        # Allow not running nucflag.
        infile_stream=lambda wc, input: (
            f"join {input.interm_bed} {input.nucflag_bed}"
            if input.nucflag_bed
            else f"cat {input.interm_bed}"
        ),
    log:
        "logs/get_complete_correct_cens/get_complete_correct_bed_{sm}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ {params.infile_stream} | \
        awk -v OFS="\\t" '{{
            status=(($8 == "good" || $8 == "") && ($4 == "false")) ? "good" : "misassembled";
            if (status != "good") {{
                next;
            }};
            split($1, arr, ":");
            print arr[1], $2, $3
        }}' ;}} > {output} 2> {log}
        """


use rule extract_and_index_fa as get_complete_correct_cens_fa with:
    input:
        # fasta with no_coords
        fa=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}_regions.renamed.reort.fa",
        ),
        # (name_no_coords, st, end)
        bed=rules.get_complete_correct_cens_bed.output,
    output:
        seq=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "seq",
            "final",
            "{sm}_cens.fa",
        ),
        idx=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "seq",
            "final",
            "{sm}_cens.fa.fai",
        ),
    log:
        "logs/get_complete_correct_cens/get_complete_correct_cens_fa_{sm}.log",
    params:
        added_cmds="",
        bed=lambda wc, input: f"""<(awk -v OFS="\\t" '{{print $1, $2-1, $3}}' {input.bed})""",


rule get_complete_correct_cens_all:
    input:
        expand(rules.get_complete_correct_cens_bed.output, sm=SAMPLE_NAMES),
        expand(rules.get_complete_correct_cens_fa.output, sm=SAMPLE_NAMES),
