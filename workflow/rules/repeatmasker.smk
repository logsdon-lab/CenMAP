
rule refmt_correct_alr_regions_for_rm:
    input:
        seq=lambda wc: expand(
            rules.extract_alr_region_sample_by_chr.output.seq,
            ort=ORIENTATION,
            chr=CHROMOSOMES,
            sm=[wc.sm],
        ),
        bed=get_correct_assembly_regions,
    output:
        seq=temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "{sm}_joined_ALR_regions.500kbp.fa",
            )
        ),
        lst=temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "{sm}_ALR_regions.500kbp.lst",
            )
        ),
    log:
        "logs/refmt_alr_regions_repeatmasker_{sm}.log",
    shell:
        """
        # Join fwd and rev alr regions fa.
        cat {input.seq} > {output.seq}
        # (sm)_(chr)_(ctg):(start)-(end)
        awk -v OFS="\\t" '{{ print "{wildcards.sm}_"$5"_"$1":"$2"-"$3 }}' {input.bed} > {output.lst} 2> {log}
        """


use rule extract_and_index_fa as extract_correct_alr_regions_rm with:
    input:
        fa=rules.refmt_correct_alr_regions_for_rm.output.seq,
        bed=rules.refmt_correct_alr_regions_for_rm.output.lst,
    output:
        seq=os.path.join(
            config["repeatmasker"]["output_dir"], "{sm}_correct_ALR_regions.500kbp.fa"
        ),
        idx=os.path.join(
            config["repeatmasker"]["output_dir"],
            "{sm}_correct_ALR_regions.500kbp.fa.fai",
        ),
    params:
        added_cmds="",
    log:
        "logs/extract_alr_regions_repeatmasker_{sm}.log",


# TODO: Doc uses modified version. May need to update ID length from 50 -> 70 so ask glennis.
# https://github.com/search?q=repo%3Armhubley%2FRepeatMasker%20%2050&type=code
rule run_repeatmasker:
    input:
        script="workflow/scripts/RepeatMasker/RepeatMasker",
        seq=rules.extract_correct_alr_regions_rm.output.seq,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "{sm}",
            "{sm}_correct_ALR_regions.500kbp.fa.out",
        ),
    conda:
        "../env/tools.yaml"
    threads: config["repeatmasker"]["threads"]
    params:
        output_dir=lambda wc, output: os.path.dirname(str(output)),
        species="human",
        engine="rmblast",
    log:
        "logs/repeatmasker_{sm}.log",
    benchmark:
        "benchmarks/repeatmasker_{sm}.tsv"
    shell:
        """
        perl {input.script} \
        -engine {params.engine} \
        -species {params.species} \
        -dir {params.output_dir} \
        -pa {threads} \
        {input.seq} &> {log}
        """


# TODO: Grp by chr


# Run repeatmasker on reference t2t-chm13 as a control.
use rule run_repeatmasker as run_repeatmasker_ref with:
    input:
        seq=config["align_asm_to_ref"]["config"]["ref"][REF_NAME],
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            REF_NAME,
            f"{REF_NAME}_correct_ALR_regions.500kbp.fa.out",
        ),
    log:
        f"logs/repeatmasker_{REF_NAME}.log",
    benchmark:
        f"benchmarks/repeatmasker_{REF_NAME}.tsv"


# Reformat repeatmasker
# |1259|28.4|7.4|5.3|GM18989_chr1_hap1-0000003:9717731-15372230|8|560|(5653940)|+|Charlie2b|DNA/hAT-Charlie|120|683|(2099)|1|
rule format_repeatmasker_output:
    input:
        expand(
            rules.run_repeatmasker.output,
            sm=SAMPLE_NAMES,
        ),
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "all_samples_correct_ALR_regions.500kbp.fa.out",
        ),
    log:
        "logs/format_repeatmasker_output_all_samples.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ cat {input} | \
        sed -e 's/:/\\t/g' -e 's/-/\\t/g' | \
        awk -v OFS="\\t" '{{print $1, $2, $3, $4, $5"-"$6":"$7"-"$8, $9, $10, $11, $12, $13,
$14, $15, $16, $17, $18}}' | \
        tail -n +4;}} > {output} 2> {log}
        """


rule format_add_control_repeatmasker_output:
    input:
        ref_rm_output=rules.run_repeatmasker_ref.output,
        sample_rm_output=rules.format_repeatmasker_output.output,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "all_samples_and_ref_correct_ALR_regions.500kbp.fa.out",
        ),
    log:
        "logs/format_add_control_repeatmasker_output.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        # Copy file and append reference repeatmasker output.
        cp {input.sample_rm_output} {output} 2> {log}
        {{ grep "chr" {input.ref_rm_output} | \
        sed -e 's/chr/chm13_chr/g' -e 's/:/\\t/g' -e 's/-/\\t/g' | \
        awk -v OFS="\\t" \
        '{{print $1, $2, $3, $4, $5":"$6"-"$7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17}}';}} >> {output} 2> {log}
        """


rule format_add_censat_annot_repeatmasker_output:
    input:
        # TODO: Where at? Should have rule for this as well?
        cen_sat_rm_output=config["repeatmasker"]["censat_annot_hor_output"],
        sample_ctrl_rm_output=rules.format_add_control_repeatmasker_output.output,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "all_correct_ALR_regions.500kbp.fa.out",
        ),
    log:
        "logs/format_add_censat_annot_repeatmasker_output.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        cp {input.sample_ctrl_rm_output} {output} 2> {log}

        {{ cat {input.cen_sat_rm_output} | \
        sed -e 's/chr/chm13_chr/g' -e 's/:/\\t/g' -e 's/-/\\t/g' | \
        awk -v OFS="\\t" '{{print $1, $2, $3,
$4, $5":"$6"-"$7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17}}';}} >> {output} 2> {log}
        """


rule reverse_complete_repeatmasker_output:
    input:
        rules.format_add_censat_annot_repeatmasker_output.output,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "all_correct_ALR_regions.500kbp.rc.fa.out",
        ),
    log:
        "logs/reverse_complete_repeatmasker_output.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        tac {input} | \
        sed -e 's/:/\\t/g' -e 's/-/\\t/g' | \
        awk -v OFS="\\t" '{{print $1, $2, $3, $4, $5"-"$6":"$7"-"$8, $8-$7-$10+1,
$8-$7-$9+1, $11, $12, $13, $14, $15, $16, $17, $18}}' | \
        sed 's/chr/rc_chr/g' | \
        grep -v "chm13" > {output}
        """
