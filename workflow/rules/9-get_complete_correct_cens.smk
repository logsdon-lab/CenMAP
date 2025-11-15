include: "common.smk"
include: "utils.smk"
include: "7-fix_cens_w_repeatmasker.smk"
include: "8-humas_annot.smk"

COMPLETE_CORRECT_CENS_LOGDIR = join(LOG_DIR, "9-get_complete_correct_cens")


if config.get("nucflag"):

    include: "8-nucflag.smk"


rule get_complete_correct_cens_bed_unfiltered:
    input:
        # (name, st, end, ctg_name, ctg_len)
        interm_bed=ancient(rules.make_complete_cens_bed.output.cen_bed),
        # (name, st, end, status)
        nucflag_bed=(
            rules.check_asm_nucflag.output.asm_status if config.get("nucflag") else []
        ),
    output:
        # BED9
        # (name, st, end, adj_name, score, ort, adj_st, adj_end, rgb)
        bed=temp(join(
            FINAL_OUTPUT_DIR,
            "bed",
            "{sm}_complete_correct_cens.unfiltered.bed",
        )),
    params:
        # Allow not running nucflag.
        infile_stream=lambda wc, input: (
            f"join <(sort -k1 {input.interm_bed}) <(sort -k1 {input.nucflag_bed})"
            if input.nucflag_bed
            else f"cat {input.interm_bed}"
        ),
    log:
        join(COMPLETE_CORRECT_CENS_LOGDIR, "get_complete_correct_bed_unfiltered_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ {params.infile_stream} | \
        sort -k1,1 -k2,2n | \
        uniq | \
        awk -v OFS="\\t" '{{
            adj_name=$1;
            adj_st=$2; adj_end=$3;
            ctg_name=$4;
            ctg_len=$5;
            ort=($1 ~ "rc-") ? "-" : "+";
            ctg_st=$2; ctg_end=$3;
            if (ort == "-") {{
                new_ctg_st=ctg_len-ctg_end + 1;
                new_ctg_end=ctg_len-ctg_st + 1;
                ctg_st=new_ctg_st;
                ctg_end=new_ctg_end;
            }}
            status=($8 == "good" || $8 == "") ? "good" : "misassembled";
            if (status != "good") {{
                next;
            }};
            print ctg_name, ctg_st, ctg_end, adj_name, 0, ort, adj_st, adj_end, "0,0,0"
        }}' ;}} | \
        sort -k1,1 -k2,2n | \
        uniq \
        > {output} 2> {log}
        """

# Could not reuse the rule from 10-format_hor_stv.smk due to circular import issues.
def humas_annot_all_outputs(wc):
    _ = [checkpoints.run_humas_annot.get(sm=sm).output for sm in SAMPLE_NAMES]
    fastas = glob.glob(join(HUMAS_CENS_SPLIT_DIR, "*.fa"))
    fnames = get_valid_fnames(fastas, filter_chrom = None)
    outputs = config.get("plot_hor_stv", {}).get("ref_stv", [])
    if config["humas_annot"]["mode"] == "srf-n-trf":
        outputs.extend(
            expand(
                rules.filter_srf_trf_annot.output,
                zip,
                sm=[fname.sm for fname in fnames],
                fname=fnames,
            )
        )
    elif config["humas_annot"]["mode"] == "sf":
        outputs.extend(
            expand(rules.format_monomer_sf_classes.output, fname=fnames, mode="sf")
        )
    else:
        return expand(rules.cens_generate_stv.output, zip, fname=fnames, chr=CHROMOSOMES if CHROMOSOMES else ["all"])

    return outputs

rule merge_all_stv_bed:
    input:
        humas_annot_all_outputs,
    output:
        merged_stv=temp(join(
            FINAL_OUTPUT_DIR,
            "bed",
            "all_stv_merged.bed",
        )),
    log:
        join(COMPLETE_CORRECT_CENS_LOGDIR, "merge_all_stv_bed.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        cat {input} \
        | cut -f1 \
        | sort \
        | uniq \
         > {output.merged_stv} 2> {log}
        """

rule get_complete_correct_cens_bed:
    input:
        cens_bed_unfiltered = rules.get_complete_correct_cens_bed_unfiltered.output,
        merged_stv_bed = rules.merge_all_stv_bed.output.merged_stv if config.get("get_complete_correct_cens", {}).get("filter_by_hor", True) else "", 
    output:
        # BED9
        # (name, st, end, adj_name, score, ort, adj_st, adj_end, rgb)
        bed=join(
            FINAL_OUTPUT_DIR,
            "bed",
            "{sm}_complete_correct_cens.bed",
        ),
    log:
        join(COMPLETE_CORRECT_CENS_LOGDIR, "get_complete_correct_bed_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        # Filter cens_bed_unfiltered by matching col1:col2-col3 against patterns
        awk -v OFS="\\t" '
        NR==FNR {{
            patterns[$1] = 1;
            next;
        }}
        {{
            key = $4 ":" ($7+1) "-" $8;
            if (key in patterns) {{
                print $0;
            }}
        }}' <(cat {input.merged_stv_bed} | grep "^{wildcards.sm}_") {input.cens_bed_unfiltered} > {output.bed}
        2> {log}
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
