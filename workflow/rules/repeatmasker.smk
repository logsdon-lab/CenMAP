
rule get_valid_regions_for_rm:
    input:
        bed=get_correct_assembly_regions,
    output:
        good_regions=os.path.join(
            config["repeatmasker"]["output_dir"],
            "{sm}_trimmed_correct_ALR_regions.500kbp.bed",
        ),
    shell:
        """
        grep "good" {input.bed} > {output.good_regions}
        """


use rule extract_and_index_fa as extract_correct_alr_regions_rm with:
    input:
        fa=rules.concat_asm.output,
        bed=rules.get_valid_regions_for_rm.output,
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


rule count_complete_cens:
    input:
        expand(rules.extract_correct_alr_regions_rm.output.idx, sm=SAMPLE_NAMES),
    output:
        all_idx=os.path.join(
            config["repeatmasker"]["output_dir"],
            "all_correct_ALR_regions.500kbp.fa.fai",
        ),
        cmp_cnts=os.path.join(
            config["repeatmasker"]["output_dir"],
            "complete_cen_counts.tsv",
        ),
    params:
        strip_dir=config["repeatmasker"]["output_dir"].replace("/", "\/") + "\/",
    shell:
        """
        rm {output.all_idx}
        cat {input} > {output.all_idx}

        wc -l {input} {output.all_idx} | \
        sed -e 's/{params.strip_dir}//g' -e 's/_/\\t/g' | \
        awk -v OFS="\\t" '{{print $2, $1, $1/46*100}}' | head -n +2 > {output}
        """


# TODO: Doc uses modified version. May need to update ID length from 50 -> 70 so ask glennis.
# https://github.com/search?q=repo%3Armhubley%2FRepeatMasker%20%2050&type=code
rule run_repeatmasker:
    input:
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
        RepeatMasker \
        -engine {params.engine} \
        -species {params.species} \
        -dir {params.output_dir} \
        -pa {threads} \
        {input.seq} &> {log}
        """


rule merge_legends_for_rm:
    input:
        lambda wc: expand(
            rules.new_cens_create_oriented_ctg_name_legend.output,
            sm=[wc.sm],
            ort=ORIENTATION,
        ),
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"], temp("{sm}_merged_legend.txt")
        ),
    shell:
        """
        cat {input} > {output}
        """


# TODO: use legend to replace contig name in
rule rename_contig_name_repeatmasker:
    input:
        rm_out=rules.run_repeatmasker.output,
        legend=rules.merge_legends_for_rm.output,
    output:
        renamed_out=os.path.join(
            config["repeatmasker"]["output_dir"],
            "{sm}",
            "{sm}_renamed_correct_ALR_regions.500kbp.fa.out",
        ),
    run:
        with (
            open(input.legend, "rt") as legend_fh,
            open(input.rm_out, "rt") as rm_out_fh,
            open(output.renamed_out, "rt") as rm_renamed_out_fh,
        ):
            rm_out = rm_out_fh.read()

            contig_name_legend = {
                dict(line.strip().split("\t")) for line in legend_fh.readlines()
            }
            for contig, new_name in rm_out.items():
                rm_out = rm_out.replace(contig, new_name)

            rm_renamed_out_fh.write(rm_out)


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
            rules.rename_contig_name_repeatmasker.output,
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
