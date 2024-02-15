
rule get_valid_regions_for_rm:
    input:
        bed=get_correct_assembly_regions,
    output:
        good_regions=os.path.join(
            config["repeatmasker"]["output_dir"],
            "{sm}_valid_ALR_regions.500kbp.bed",
        ),
    params:
        assembly_filter="good",
    shell:
        """
        grep "{params.assembly_filter}" {input.bed} > {output.good_regions}
        """


rule merge_ort_alr_regions_for_rm:
    input:
        lambda wc: expand(
            rules.extract_new_oriented_cens_regions.output.seq,
            sm=[wc.sm],
            ort=ORIENTATION,
        ),
    output:
        temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "{sm}_merged_ALR_regions.500kbp.fa",
            )
        ),
    shell:
        """
        cat {input} > {output}
        """


use rule extract_and_index_fa as extract_correct_alr_regions_rm with:
    input:
        fa=rules.merge_ort_alr_regions_for_rm.output,
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
        # dir to remove from output strings.
        strip_dir=config["repeatmasker"]["output_dir"].replace("/", "\/") + "\/",
        num_chrs=46,
    log:
        "logs/count_complete_cens.log",
    shell:
        """
        rm -f {output.all_idx}
        cat {input} > {output.all_idx}

        {{ wc -l {input} {output.all_idx} | \
        sed -e 's/{params.strip_dir}//g' -e 's/_/\\t/g' | \
        awk -v OFS="\\t" '{{print $2, $1, $1/{params.num_chrs}*100}}' | \
        head -n +2;}} > {output.cmp_cnts} 2> {log}
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


# TODO: rules.format_repeatmasker_output can be joined here probably.
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
        import csv

        with (
            open(str(input.legend), "rt") as legend_fh,
            open(str(input.rm_out), "rt") as rm_out_fh,
            open(str(output.renamed_out), "wt") as rm_renamed_out_fh,
        ):
            contig_name_legend = dict(
                tuple(line.strip().split("\t")) for line in legend_fh.readlines()
            )
            headers = [next(rm_out_fh) for _ in range(3)]
            rm_renamed_out_writer = csv.writer(rm_renamed_out_fh, delimiter="\t")
            # Write headers to file handle without csv writer to not screw up formatting.
            rm_renamed_out_fh.writelines(headers)

            for rm_row in rm_out_fh.readlines():
                rm_row = rm_row.strip().split()
                # From: {ctg}:{start}-{end}
                # To: {sm}_{chr}_{ctg}
                query_seq_name, coords = rm_row[4].split(":")
                new_query_seq_name = contig_name_legend[query_seq_name]
                rm_row[4] = new_query_seq_name
                rm_renamed_out_writer.writerow(rm_row)


# Run repeatmasker on reference t2t-chm13 as a control.
# TODO: Do I need to format this like below?
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
        ref_rm_output=(
            config["repeatmasker"]["ref_repeatmasker_output"]
            if config["repeatmasker"].get("ref_repeatmasker_output")
            else rules.run_repeatmasker_ref.output
        ),
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
        {{ tac {input} | \
        sed -e 's/:/\\t/g' -e 's/-/\\t/g' | \
        awk -v OFS="\\t" '{{print $1, $2, $3, $4, $5"-"$6":"$7"-"$8, $8-$7-$10+1,
$8-$7-$9+1, $11, $12, $13, $14, $15, $16, $17, $18}}' | \
        sed 's/chr/rc_chr/g' | \
        grep -v "chm13";}} > {output} 2> {log}
        """


rule plot_cens_from_rm_by_chr:
    input:
        script="workflow/scripts/repeatStructure_onlyRM.R",
        rm_out=rules.format_add_censat_annot_repeatmasker_output.output,
    output:
        rm_out_by_chr=os.path.join(
            config["repeatmasker"]["output_dir"], "{chr}_cens.fa.out"
        ),
        repeat_plot_by_chr=os.path.join(
            config["repeatmasker"]["output_dir"], "hgsvc3_{chr}_cens.additional.pdf"
        ),
    log:
        "logs/plot_{chr}_cens_from_rm.log",
    conda:
        "../env/r.yaml"
    shell:
        """
        grep "{wildcards.chr}_" {input.rm_out} > {output.rm_out_by_chr}
        grep "{wildcards.chr}:" {input.rm_out} >> {output.rm_out_by_chr}
        Rscript {input.script} {output.rm_out_by_chr} {output.repeat_plot_by_chr} 2> {log}
        """
