
include: "common.smk"
include: "utils.smk"


wildcard_constraints:
    chr="|".join(CHROMOSOMES),


rule get_valid_regions_for_rm:
    input:
        bed=(
            os.path.join(
                config["nucflag"]["output_dir"],
                "{sm}_cen_status.bed",
            )
            if "nucflag" in config
            else os.path.join(
                config["new_cens"]["output_dir"], "bed", "{sm}_ALR_regions.bed"
            )
        ),
    output:
        good_regions=os.path.join(
            config["repeatmasker"]["output_dir"],
            "bed",
            "{sm}_valid_ALR_regions.bed",
        ),
    params:
        omit_nucflag="nucflag" not in config,
        assembly_filter="good",
    shell:
        """
        awk -v OFS="\\t" '{{
            if ("{params.omit_nucflag}" == "True") {{
                $4="good"
            }}
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
            "{sm}_regions.renamed.reort.fa",
        ),
        bed=rules.get_valid_regions_for_rm.output,
    output:
        seq=temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "seq",
                "{sm}_correct_ALR_regions.fa",
            )
        ),
        idx=temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "seq",
                "{sm}_correct_ALR_regions.fa.fai",
            )
        ),
    log:
        "logs/repeatmasker/extract_alr_regions_repeatmasker_{sm}.log",


# RepeatMasker has a limit of 50 characters for sequence names.
# While I was able to create a fork of RepeatMasker that allowed for longer sequence names, I still ran into issues.
# Bumping to v4.1.5 and also adding the increased limit fixed the issue but resulted in a 5-10x increase in runtime.
# So I downgraded the repeatmasker70 image to v4.1.0.
rule rename_for_repeatmasker:
    input:
        fa=rules.extract_correct_alr_regions_rm.output.seq,
    output:
        renamed_fa=temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "seq",
                "{sm}_correct_ALR_regions.renamed.fa",
            )
        ),
        renamed_fa_idx=temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "seq",
                "{sm}_correct_ALR_regions.renamed.fa.fai",
            )
        ),
    params:
        prefix="seq",
    conda:
        "../env/tools.yaml"
    log:
        "logs/repeatmasker/rename_for_repeatmasker_{sm}.log",
    shell:
        """
        seqtk rename {input.fa} {params.prefix} > {output.renamed_fa} 2> {log}
        if [ -s {output.renamed_fa} ]; then
            samtools faidx {output.renamed_fa} 2>> {log}
        else
            touch {output.renamed_fa_idx}
        fi
        """


rule run_repeatmasker:
    input:
        seq=rules.rename_for_repeatmasker.output.renamed_fa,
    output:
        temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "repeats",
                "{sm}",
                "{sm}_correct_ALR_regions.renamed.fa.out",
            )
        ),
    threads: config["repeatmasker"]["threads"]
    params:
        output_dir=lambda wc, output: os.path.dirname(str(output)),
        species="human",
        engine="rmblast",
    conda:
        "../env/repeatmasker.yaml"
    log:
        "logs/repeatmasker/repeatmasker_{sm}.log",
    # Retry in case of .RepeatMaskerCache failure.
    retries: 2
    benchmark:
        "benchmarks/repeatmasker/repeatmasker_{sm}.tsv"
    shell:
        """
        if [ -s {input.seq} ]; then
            RepeatMasker \
            -engine {params.engine} \
            -species {params.species} \
            -dir {params.output_dir} \
            -pa {threads} \
            {input.seq} &> {log}
        else
            touch {output}
        fi
        """


# Rename repeatmasker output to match the original sequence names.
rule reformat_repeatmasker_output:
    input:
        script="workflow/scripts/reformat_rm.py",
        rm_out=rules.run_repeatmasker.output,
        original_fai=rules.extract_correct_alr_regions_rm.output.idx,
        renamed_fai=rules.rename_for_repeatmasker.output.renamed_fa_idx,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "{sm}",
            "{sm}_correct_ALR_regions.fa.reformatted.out",
        ),
    log:
        "logs/repeatmasker/reformat_repeatmasker_output_{sm}.log",
    conda:
        "../env/py.yaml"
    shell:
        """
        python {input.script} -i {input.rm_out} -of {input.original_fai} -rf {input.renamed_fai} > {output} 2> {log}
        """


# Merge repeatmasker and convert sep to tab.
# |1259|28.4|7.4|5.3|GM18989_chr1_hap1-0000003:9717731-15372230|8|560|(5653940)|+|Charlie2b|DNA/hAT-Charlie|120|683|(2099)|1|
rule merge_repeatmasker_output:
    input:
        expand(rules.reformat_repeatmasker_output.output, sm=SAMPLE_NAMES),
    output:
        temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "repeats",
                "all",
                "all_samples_correct_ALR_regions.fa.out",
            )
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
        temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "repeats",
                "ref",
                "ref_ALR_regions.fa.out",
            )
        ),


rule format_add_control_repeatmasker_output:
    input:
        ref_rm_output=rules.merge_control_repeatmasker_output.output,
        sample_rm_output=rules.merge_repeatmasker_output.output,
    output:
        temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "repeats",
                "all",
                "all_samples_and_ref_correct_ALR_regions.fa.out",
            )
        ),
    log:
        "logs/repeatmasker/format_add_control_repeatmasker_output.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        # Copy file and append reference repeatmasker output.
        cp {input.sample_rm_output} {output} 2> {log}
        grep "chr" {input.ref_rm_output} >> {output} 2> {log}
        """


rule extract_rm_out_by_chr:
    input:
        rm_out=rules.format_add_control_repeatmasker_output.output,
    output:
        rm_out_by_chr=temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "repeats",
                "all",
                "{chr}_cens.fa.out",
            )
        ),
    log:
        "logs/repeatmasker/extract_{chr}_cens_from_rm.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ grep "{wildcards.chr}[_:]" {input.rm_out} || true; }}> {output.rm_out_by_chr} 2> {log}
        """


include: "fix_cens_w_repeatmasker.smk"


rule repeatmasker_only:
    input:
        expand(rules.get_complete_correct_cens_bed.output, sm=SAMPLE_NAMES),
        expand(rules.fix_ort_asm_final.output, sm=SAMPLE_NAMES),
        expand(rules.extract_sm_complete_correct_cens.output, sm=SAMPLE_NAMES),
        expand(rules.merge_all_complete_correct_cens_fa.output, sm=SAMPLE_NAMES),
        expand(rules.plot_cens_from_rm_by_chr.output, chr=CHROMOSOMES),
        expand(rules.plot_cens_from_original_rm_by_chr.output, chr=CHROMOSOMES),
