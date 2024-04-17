
include: "common.smk"
include: "utils.smk"


wildcard_constraints:
    chr="|".join(CHROMOSOMES),


rule get_valid_regions_for_rm:
    input:
        bed=os.path.join(
            config["nuc_freq"]["output_dir"],
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
        "logs/extract_alr_regions_repeatmasker_{sm}.log",


rule rc_correct_alr_regions_rm:
    input:
        fa=rules.extract_correct_alr_regions_rm.output.seq,
    output:
        rc_seq=os.path.join(
            config["repeatmasker"]["output_dir"],
            "seq",
            "{sm}_correct_ALR_regions.500kbp.rc.fa",
        ),
    conda:
        "../env/tools.yaml"
    log:
        "logs/rc_alr_regions_repeatmasker_{sm}.log",
    shell:
        "seqtk seq -r {input.fa} > {output.rc_seq} 2> {log}"


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
        "logs/format_add_control_repeatmasker_output.log",
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
        "logs/reverse_complete_repeatmasker_output.log",
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
        "logs/extract_{chr}_cens_from_rm.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ grep "{wildcards.chr}_" {input.rm_out} || true; }}> {output.rm_out_by_chr} 2> {log}
        {{ grep "{wildcards.chr}:" {input.rm_out} || true; }} >> {output.rm_out_by_chr} 2>> {log}
        """


rule check_cens_status:
    input:
        rm_out=rules.extract_rm_out_by_chr.output.rm_out_by_chr,
        rm_ref=rules.merge_control_repeatmasker_output.output,
    output:
        cens_status=os.path.join(
            config["repeatmasker"]["output_dir"], "status", "{chr}_cens_status.tsv"
        ),
    params:
        edge_len=100_000,
        edge_perc_alr_thr=0.7,
        dst_perc_thr=0.3,
        # Edge-case for chrs whose repeats are small and broken up.
        max_alr_len_thr=lambda wc: 0 if wc.chr in ["chrY", "chr11", "chr8"] else 200_000,
    log:
        "logs/check_cens_status_{chr}.log",
    conda:
        "../env/cen_stats.yaml"
    shell:
        """
        cen-stats \
        -i {input.rm_out} \
        -r {input.rm_ref} \
        -o {output} \
        --edge_len {params.edge_len} \
        --edge_perc_alr_thr {params.edge_perc_alr_thr} \
        --dst_perc_thr {params.dst_perc_thr} \
        --max_alr_len_thr {params.max_alr_len_thr} 2> {log}
        """


rule create_correct_oriented_cens:
    input:
        rm_chr_out=rules.extract_rm_out_by_chr.output.rm_out_by_chr,
        rm_rev_out=rules.reverse_complete_repeatmasker_output.output,
        cens_correction_list=rules.check_cens_status.output.cens_status,
    output:
        reoriented_rm_out=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "reoriented_{chr}_cens.fa.out",
        ),
    log:
        "logs/create_correct_oriented_{chr}_cens_list.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        contigs_to_reverse=$(awk '{{ if ($3=="rev") {{ print $1 }} }}' {input.cens_correction_list} | sort | uniq)
        joined_contigs_to_reverse=$(echo "${{contigs_to_reverse[*]}}" | tr '\\n' '|' | sed 's/.$//')
        joined_contigs_to_reverse_rc=$(echo "${{joined_contigs_to_reverse[@]}}" | sed 's/chr/rc-chr/g' )

        # If nothing to reverse, just copy file.
        if [ -z $joined_contigs_to_reverse ]; then
            cp {input.rm_chr_out} {output.reoriented_rm_out}
        else
            # non-matching so everything correctly oriented
            grep -vE "$joined_contigs_to_reverse" {input.rm_chr_out} > {output.reoriented_rm_out}
            grep -E "$joined_contigs_to_reverse_rc" {input.rm_rev_out} >> {output.reoriented_rm_out}
        fi
        """


rule merge_corrections_list:
    input:
        expand(rules.check_cens_status.output.cens_status, chr=CHROMOSOMES),
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"], "status", "all_cen_corrections.tsv"
        ),
    shell:
        """
        cat {input} > {output}
        """


rule fix_incorrect_merged_legend:
    input:
        script="workflow/scripts/fix_incorrect_merged_legend.py",
        cens_correction_list=rules.merge_corrections_list.output,
        merged_legend=os.path.join(
            config["concat_asm"]["output_dir"], "{sm}", "{sm}_legend.txt"
        ),
    output:
        corrected_legend=os.path.join(
            config["repeatmasker"]["output_dir"],
            "status",
            "{sm}_corrected_merged_legend.txt",
        ),
    conda:
        "../env/py.yaml"
    log:
        "logs/fix_incorrect_merged_legend_{sm}.log",
    shell:
        """
        python {input.script} \
        -ic {input.cens_correction_list} \
        -l {input.merged_legend} > {output}
        """


rule fix_incorrect_mapped_cens:
    input:
        script="workflow/scripts/fix_incorrect_mapped_cens.py",
        cens_correction_list=rules.merge_corrections_list.output,
        reoriented_rm_out=expand(
            rules.create_correct_oriented_cens.output, chr=CHROMOSOMES
        ),
    output:
        corrected_cens_list=os.path.join(
            config["repeatmasker"]["output_dir"], "status", "all_corrected_cens.list"
        ),
        corrected_rm_out=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "all_corrected_cens.fa.out",
        ),
    conda:
        "../env/py.yaml"
    log:
        "logs/fix_incorrect_mapped_cens.log",
    shell:
        """
        python {input.script} \
        -ic {input.cens_correction_list} \
        -ir {input.reoriented_rm_out} > {output.corrected_rm_out} 2> {log}

        cut -f 5 {output.corrected_rm_out} | sort | uniq > {output.corrected_cens_list}
        """


rule count_complete_cens:
    input:
        rules.fix_incorrect_mapped_cens.output.corrected_cens_list,
    output:
        cmp_cnts=os.path.join(
            config["repeatmasker"]["output_dir"],
            "status",
            "complete_cen_counts.tsv",
        ),
    params:
        num_chrs=46,
    run:
        from collections import Counter

        with (
            open(str(input)) as cens_list_fh,
            open(str(output), "wt") as cen_counts_fh,
        ):
            counts = Counter()
            for l in cens_list_fh.readlines():
                if l.startswith("chm13"):
                    continue
                sample, _, _ = l.rsplit("_", maxsplit=2)
                counts[sample] += 1

            for sm, count in counts.items():
                prop = (count / params.num_chrs) * 100
                cen_counts_fh.write(f"{sm}\t{count}\t{prop}\n")

            total_count = sum(counts.values())
            total_possible_count = len(counts) * params.num_chrs
            total_prop = (total_count / total_possible_count) * 100
            cen_counts_fh.write(f"all\t{total_count}\t{total_prop}\n")


rule split_corrected_rm_output:
    input:
        corrected_cens_list=rules.fix_incorrect_mapped_cens.output.corrected_cens_list,
        corrected_rm_out=rules.fix_incorrect_mapped_cens.output.corrected_rm_out,
    output:
        corrected_rm_out=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "corrected_{chr}_cens.fa.out",
        ),
        corrected_cens_list=os.path.join(
            config["repeatmasker"]["output_dir"], "status", "corrected_{chr}_cens.list"
        ),
    log:
        "logs/split_corrected_{chr}_rm_output.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ grep "{wildcards.chr}[_:]" {input.corrected_rm_out} || true; }} > {output.corrected_rm_out}
        {{ grep "{wildcards.chr}[_:]" {input.corrected_cens_list} || true; }} > {output.corrected_cens_list}
        """


rule plot_cens_from_rm_by_chr:
    input:
        script="workflow/scripts/repeatStructure_onlyRM.R",
        rm_out=rules.split_corrected_rm_output.output.corrected_rm_out,
    output:
        repeat_plot_by_chr=os.path.join(
            config["repeatmasker"]["output_dir"],
            "plot",
            "{chr}_cens.corrected.pdf",
        ),
    log:
        "logs/plot_{chr}_cens_from_rm.log",
    conda:
        "../env/r.yaml"
    shell:
        """
        Rscript {input.script} {input.rm_out} {output.repeat_plot_by_chr} 2> {log}
        """


use rule plot_cens_from_rm_by_chr as plot_og_cens_from_rm_by_chr with:
    input:
        script="workflow/scripts/repeatStructure_onlyRM.R",
        rm_out=rules.extract_rm_out_by_chr.output.rm_out_by_chr,
    output:
        repeat_plot_by_chr=os.path.join(
            config["repeatmasker"]["output_dir"],
            "plot",
            "{chr}_cens.original.pdf",
        ),
    log:
        "logs/plot_{chr}_cens_from_rm_og.log",


rule repeatmasker_only:
    input:
        expand(rules.extract_correct_alr_regions_rm.output, sm=SAMPLE_NAMES),
        expand(rules.rc_correct_alr_regions_rm.output, sm=SAMPLE_NAMES),
        expand(rules.run_repeatmasker.output, sm=SAMPLE_NAMES),
        expand(rules.remove_repeatmasker_header.output, sm=SAMPLE_NAMES),
        rules.merge_repeatmasker_output.output,
        rules.merge_control_repeatmasker_output.output,
        rules.format_add_control_repeatmasker_output.output,
        rules.reverse_complete_repeatmasker_output.output,
        expand(rules.extract_rm_out_by_chr.output, chr=CHROMOSOMES),
        expand(rules.check_cens_status.output, chr=CHROMOSOMES),
        rules.merge_corrections_list.output,
        expand(rules.fix_incorrect_merged_legend.output, sm=SAMPLE_NAMES),
        expand(rules.fix_incorrect_mapped_cens.output, chr=CHROMOSOMES),
        rules.count_complete_cens.output,
        expand(rules.create_correct_oriented_cens.output, chr=CHROMOSOMES),
        expand(rules.fix_incorrect_mapped_cens.output, chr=CHROMOSOMES),
        expand(rules.split_corrected_rm_output.output, chr=CHROMOSOMES),
        expand(rules.plot_cens_from_rm_by_chr.output, chr=CHROMOSOMES),
        expand(rules.plot_og_cens_from_rm_by_chr.output, chr=CHROMOSOMES),
