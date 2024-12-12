
rule get_complete_correct_cens_bed:
    input:
        # (name, st, end, is_partial, ctg_name, ctg_len)
        interm_bed=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            "interm",
            "{sm}_complete_cens.bed",
        ),
        # (name, st, end, status)
        nucflag_bed=(
            os.path.join(
                config["nucflag"]["output_dir"],
                "{sm}_cen_status.bed",
            )
            if config.get("nucflag")
            else []
        ),
    output:
        # BED9
        # (name, st, end, adj_name, score, ort, adj_st, adj_end, rgb)
        bed=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            "final",
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
        "logs/get_complete_correct_cens/get_complete_correct_bed_{sm}.log",
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
                new_ctg_st=ctg_len-ctg_end; new_ctg_end=ctg_len-ctg_st;
                ctg_st=new_ctg_st; ctg_end=new_ctg_end;
            }}
            status=(($9 == "good" || $9 == "") && ($4 == "false")) ? "good" : "misassembled";
            if (status != "good") {{
                next;
            }};
            print ctg_name, ctg_st, ctg_end, adj_name, 0, ort, adj_st, adj_end, "0,0,0"
        }}' ;}} > {output} 2> {log}
        """


use rule extract_and_index_fa as get_complete_correct_cens_fa with:
    input:
        # fasta with no_coords
        fa=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}-asm-renamed-reort.fa",
        ),
        # (name, st, end)
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
        bed=lambda wc, input: f"""<(awk -v OFS="\\t" '{{print $4, $7-1, $8}}' {input.bed})""",


rule get_complete_correct_cens_all:
    input:
        expand(rules.get_complete_correct_cens_bed.output, sm=SAMPLE_NAMES),
        expand(rules.get_complete_correct_cens_fa.output, sm=SAMPLE_NAMES),
