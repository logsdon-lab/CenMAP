# Included workflow. Cannot be run separately outside of `repeatmasker.smk`.


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
