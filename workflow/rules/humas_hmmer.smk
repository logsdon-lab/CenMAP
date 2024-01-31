# This is a port of the following:
# * https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/main/hmmer-run.sh
# * https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/main/hmmer-run_SF.sh


# https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/dc5e958dfc3820fb8cc21a6ccc5cd8e5d6ef1052/hmmer-run.sh#L31
HMM_PROFILE_NAME, _ = os.path.splitext(os.path.basename(config["humas_hmmer"]["model"]))
# https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/dc5e958dfc3820fb8cc21a6ccc5cd8e5d6ef1052/hmmer-run.sh#L36
CHR_NAME_PATTERN = os.path.join(config["humas_hmmer"]["input_dir"], "{chr}.fa")
CHR_NAME = glob_wildcards(CHR_NAME_PATTERN)

# Based on mode, change output file and commands.
# https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/dc5e958dfc3820fb8cc21a6ccc5cd8e5d6ef1052/hmmer-run.sh#L62
# https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/dc5e958dfc3820fb8cc21a6ccc5cd8e5d6ef1052/hmmer-run_SF.sh#L62
# https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/dc5e958dfc3820fb8cc21a6ccc5cd8e5d6ef1052/hmmer-run.sh#L66
# https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/dc5e958dfc3820fb8cc21a6ccc5cd8e5d6ef1052/hmmer-run_SF.sh#L66
HUMAS_HMMER_MODE = config["humas_hmmer"]["mode"]
if HUMAS_HMMER_MODE == "AS-HOR":
    overlap_output_pattern = "AS-HOR+SF-vs-{chr}.bed"
    final_output_pattern = "AS-HOR-vs-{chr}.bed"
    # FEDOR: AS-HOR only (skip SF monomers)
    final_output_cmd = "'{{ if (length($4)==2) {{next}} print}}'"
elif HUMAS_HMMER_MODE == "SF":
    overlap_output_pattern = "AS-SF-vs-{chr}.bed"
    final_output_pattern = "AS-strand-vs-{chr}.bed"
    # FEDOR: AS-strand annotation. "+" is blue, "-" is red
    final_output_cmd = """
    -F $'\\t' 'BEGIN {{OFS = FS}} {{if ($6=="+") {{$9="0,0,255"}}; if ($6=="-") {{$9="255,0,0"}} print $0}}'
    """
else:
    raise ValueError(f"Invalid HumAS-HMMER mode. ({config['humas_hmmer']['mode']})")


rule humas_hmmer_all:
    input:
        expand(overlap_output_pattern, chr=CHR_NAME.chr),
        expand(final_output_pattern, chr=CHR_NAME.chr),


# https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/dc5e958dfc3820fb8cc21a6ccc5cd8e5d6ef1052/hmmer-run.sh#L42
rule humas_hmmer_analysis:
    input:
        chr_fa=CHR_NAME_PATTERN,
        hmm=config["humas_hmmer"]["model"],
    output:
        # https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/dc5e958dfc3820fb8cc21a6ccc5cd8e5d6ef1052/hmmer-run.sh#L49
        # Ugh. f-string / numbers will screw up regex for wildcard parsing.
        temp("nhmmer-" + HMM_PROFILE_NAME + "-vs-{chr}-tbl.out"),
    threads: 30
    params:
        # https://manpages.ubuntu.com/manpages/focal/man1/nhmmer.1.html
        no_line_limit="--notextw",
        no_alignment="--noali",
        toss_human_readable_output="-o /dev/null",
    conda:
        "env/tools.yaml"
    log:
        "logs/humas_hmmer_analysis_{chr}.log",
    benchmark:
        "benchmarks/humas_hmmer_analysis_{chr}.tsv"
    shell:
        """
        nhmmer \
        --cpu {threads} \
        {params.no_line_limit} \
        {params.no_alignment} \
        --tblout {output} \
        {params.toss_human_readable_output} \
        {input.hmm} {input.chr_fa} &> {log}
        """


# https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/dc5e958dfc3820fb8cc21a6ccc5cd8e5d6ef1052/hmmer-run.sh#L46
rule reformat_hmm_tbl_to_bed_w_thr:
    input:
        script="workflow/scripts/hmmertblout2bed.awk",
        hmm_tbl=rules.humas_hmmer_analysis.output,
    output:
        regions=temp(HMM_PROFILE_NAME + "-vs-{chr}-tbl.bed"),
    params:
        threshold_score=0.7,
    conda:
        "env/tools.yaml"
    log:
        "logs/reformat_hmm_tbl_to_bed_w_thr_{chr}.log",
    shell:
        """
        awk \
        -v th={params.threshold_score} \
        -f {input.script} \
        {input.hmm_tbl} > {output} 2> {log}
        """


# https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/dc5e958dfc3820fb8cc21a6ccc5cd8e5d6ef1052/hmmer-run.sh#L53C4-L53C73
rule sort_hmm_res_by_name_coord:
    input:
        rules.reformat_hmm_tbl_to_bed_w_thr.output,
    output:
        temp("_nhmmer-t0-{chr}.bed"),
    conda:
        "env/tools.yaml"
    log:
        "logs/sort_hmm_res_by_name_coord_{chr}.log",
    shell:
        """
        sort -k 1.4,1 -k 2,2n {input} > {output} 2> {log}
        """


# https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/dc5e958dfc3820fb8cc21a6ccc5cd8e5d6ef1052/hmmer-run.sh#L56
rule filter_hmm_res_by_score:
    input:
        rules.sort_hmm_res_by_name_coord.output,
    output:
        temp("_nhmmer-t1-{chr}.bed"),
    # https://bedops.readthedocs.io/en/latest/content/reference/statistics/bedmap.html#usage
    params:
        # Operation
        op_opt="--max-element",
        # Overlap
        overlap_opt="--fraction-either 0.1",
    conda:
        "../env/tools.yaml"
    log:
        "logs/filter_hmm_res_by_score_{chr}.log",
    shell:
        """
        bedmap {params.op_opt} {params.overlap_opt} \
        {input} > {output} 2> {log}
        """


# https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/dc5e958dfc3820fb8cc21a6ccc5cd8e5d6ef1052/hmmer-run.sh#L60
# FEDOR: don't skip SF monomers
rule filter_hmm_res_uniq_element:
    input:
        rules.filter_hmm_res_by_score.output,
    output:
        # Original workflow overwrites sort_hmm_tbl_bed_by_name_coord output.
        # I give unique name and marked as temp.
        temp("_nhmmer-t3-{chr}.bed"),
    conda:
        "env/tools.yaml"
    log:
        "logs/filter_hmm_res_uniq_element_{chr}.log",
    shell:
        """
        awk "{{if(!(\$0 in a)){{a[\$0]; print}}}}" {input} > {output} 2> {log}
        """


rule filter_hmm_res_overlaps:
    input:
        script="workflow/scripts/overlap_filter.py",
        regions=rules.filter_hmm_res_uniq_element.output,
    output:
        overlap_output_pattern,
    conda:
        "env/py.yaml"
    log:
        "logs/filter_hmm_res_overlaps_{chr}.log",
    shell:
        """
        python {input.script} {input.regions} > {output} 2> {log}
        """


rule filter_final_hmm_res_by_mode:
    input:
        rules.filter_hmm_res_overlaps.output,
    output:
        final_output_pattern,
    params:
        cmd=lambda wc: final_output_cmd,
    conda:
        "env/tools.yaml"
    log:
        "logs/filter_final_hmm_res_by_mode_{chr}.log",
    shell:
        """
        awk {params.cmd} {input} > {output} 2> {log}
        """
