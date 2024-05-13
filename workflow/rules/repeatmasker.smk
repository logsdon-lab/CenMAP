
include: "common.smk"
include: "utils.smk"


wildcard_constraints:
    chr="|".join(CHROMOSOMES),


rule get_valid_regions_for_rm:
    input:
        bed=os.path.join(
            config["nucflag"]["output_dir"],
            "{sm}_cen_status.bed",
        ),
    output:
        good_regions=os.path.join(
            config["repeatmasker"]["output_dir"],
            "bed",
            "{sm}_valid_ALR_regions.500kbp.bed",
        ),
    params:
        assembly_filter="good",
    shell:
        """
        awk -v OFS="\\t" '{{
            if ($4 == "{params.assembly_filter}") {{
                print $1, $2, $3, $3-$2, $4
            }}
        }}' {input.bed} > {output.good_regions}
        """


use rule extract_and_index_fa as extract_correct_alr_regions_rm with:
    input:
        fa=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}",
            "{sm}_regions.renamed.fa",
        ),
        bed=rules.get_valid_regions_for_rm.output,
    output:
        seq=os.path.join(
            config["repeatmasker"]["output_dir"],
            "seq",
            "{sm}_correct_ALR_regions.500kbp.fa",
        ),
        idx=os.path.join(
            config["repeatmasker"]["output_dir"],
            "seq",
            "{sm}_correct_ALR_regions.500kbp.fa.fai",
        ),
    params:
        added_cmds="",
    log:
        "logs/repeatmasker/extract_alr_regions_repeatmasker_{sm}.log",


# rule rc_correct_alr_regions_rm:
#     input:
#         fa=rules.extract_correct_alr_regions_rm.output.seq,
#     output:
#         rc_seq=os.path.join(
#             config["repeatmasker"]["output_dir"],
#             "seq",
#             "{sm}_correct_ALR_regions.500kbp.rc.fa",
#         ),
#     conda:
#         "../env/tools.yaml"
#     log:
#         "logs/rc_alr_regions_repeatmasker_{sm}.log",
#     shell:
#         "seqtk seq -r {input.fa} > {output.rc_seq} 2> {log}"


rule run_repeatmasker:
    input:
        seq=rules.extract_correct_alr_regions_rm.output.seq,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "{sm}",
            "{sm}_correct_ALR_regions.500kbp.fa.out",
        ),
    threads: config["repeatmasker"]["threads"]
    params:
        output_dir=lambda wc, output: os.path.dirname(str(output)),
        species="human",
        engine="rmblast",
    singularity:
        "docker://logsdonlab/repeatmasker70:latest"
    log:
        "logs/repeatmasker/repeatmasker_{sm}.log",
    benchmark:
        "benchmarks/repeatmasker/repeatmasker_{sm}.tsv"
    shell:
        """
        RepeatMasker \
        -engine {params.engine} \
        -species {params.species} \
        -dir {params.output_dir} \
        -pa {threads} \
        {input.seq} &> {log}
        """


rule remove_repeatmasker_header:
    input:
        rules.run_repeatmasker.output,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "{sm}",
            "{sm}_correct_ALR_regions.500kbp.fa.noheader.out",
        ),
    shell:
        """
        tail -n +4 {input} > {output}
        """


# Merge repeatmasker and convert sep to tab.
# |1259|28.4|7.4|5.3|GM18989_chr1_hap1-0000003:9717731-15372230|8|560|(5653940)|+|Charlie2b|DNA/hAT-Charlie|120|683|(2099)|1|
rule merge_repeatmasker_output:
    input:
        expand(rules.remove_repeatmasker_header.output, sm=SAMPLE_NAMES),
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "all_samples_correct_ALR_regions.500kbp.fa.out",
        ),
    shell:
        """
        awk -v OFS="\\t" '{{$1=$1; print}}' {input} > {output}
        """


# Format repeatmasker reference output.
use rule merge_repeatmasker_output as merge_control_repeatmasker_output with:
    input:
        # Contains header. Should be first.
        config["repeatmasker"]["ref_repeatmasker_chrY_output"],
        config["repeatmasker"]["ref_repeatmasker_output"],
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "ref",
            "ref_ALR_regions.fa.out",
        ),


rule format_add_control_repeatmasker_output:
    input:
        ref_rm_output=rules.merge_control_repeatmasker_output.output,
        sample_rm_output=rules.merge_repeatmasker_output.output,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "all_samples_and_ref_correct_ALR_regions.500kbp.fa.out",
        ),
    log:
        "logs/repeatmasker/format_add_control_repeatmasker_output.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        # Copy file and append reference repeatmasker output.
        cp {input.sample_rm_output} {output} 2> {log}
        {{ grep "chr" {input.ref_rm_output} | sed -e 's/chr/chm13_chr/g';}} >> {output} 2> {log}
        """


rule reverse_complete_repeatmasker_output:
    input:
        rules.format_add_control_repeatmasker_output.output,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "all_samples_and_ref_correct_ALR_regions.500kbp.rc.fa.out",
        ),
    log:
        "logs/repeatmasker/reverse_complete_repeatmasker_output.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ tac {input} | awk -v OFS="\\t" '{{
            match($5, ":(.+)-", starts);
            match($5, ".*-(.+)$", ends);
            if ($5 ~ !/chm13/) {{
                new_start=ends[1]-starts[1]-$7+1;
                new_end=ends[1]-starts[1]-$6+1;
                $6=new_start;
                $7=new_end;
            }}
            print
        }}' | sed 's/chr/rc-chr/g' | grep -v "chm13";}} > {output} 2> {log}
        """


rule extract_rm_out_by_chr:
    input:
        rm_out=rules.format_add_control_repeatmasker_output.output,
    output:
        rm_out_by_chr=os.path.join(
            config["repeatmasker"]["output_dir"], "repeats", "all", "{chr}_cens.fa.out"
        ),
    log:
        "logs/repeatmasker/extract_{chr}_cens_from_rm.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ grep "{wildcards.chr}_" {input.rm_out} || true; }}> {output.rm_out_by_chr} 2> {log}
        {{ grep "{wildcards.chr}:" {input.rm_out} || true; }} >> {output.rm_out_by_chr} 2>> {log}
        """


include: "fix_cens_w_repeatmasker.smk"


rule repeatmasker_only:
    input:
        rules.merge_corrections_list.output,
        rules.count_complete_cens.output,
        expand(rules.plot_cens_from_rm_by_chr.output, chr=CHROMOSOMES),
        expand(rules.plot_og_cens_from_rm_by_chr.output, chr=CHROMOSOMES),
