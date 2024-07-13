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
        edge_len=lambda wc: 250_000 if wc.chr in ["chr14"] else 500_000,
        edge_perc_alr_thr=0.95,
        dst_perc_thr=0.3,
        # Edge-case for chrs whose repeats are small and broken up.
        max_alr_len_thr=lambda wc: 0 if wc.chr in ["chrY", "chr11", "chr8"] else 50_000,
        # Only allow mapping changes to 13 and 21 if chr13 or chr21.
        restrict_13_21="--restrict_13_21",
    log:
        "logs/fix_cens_w_repeatmasker/check_cens_status_{chr}.log",
    conda:
        "../env/cen_stats.yaml"
    shell:
        """
        censtats status \
        -i {input.rm_out} \
        -r {input.rm_ref} \
        -o {output} \
        --edge_len {params.edge_len} \
        --edge_perc_alr_thr {params.edge_perc_alr_thr} \
        --dst_perc_thr {params.dst_perc_thr} \
        --max_alr_len_thr {params.max_alr_len_thr} \
        {params.restrict_13_21} 2> {log}
        """


# TODO: Create second bed with rc so can extract from original assembly.
rule get_complete_correct_cens_bed:
    input:
        beds=expand(rules.get_valid_regions_for_rm.output, sm=SAMPLE_NAMES),
        statuses=expand(rules.check_cens_status.output.cens_status, chr=CHROMOSOMES),
    output:
        complete_correct_cens_bed=os.path.join(
            config["repeatmasker"]["output_dir"],
            "status",
            "all_complete_correct_cens.bed",
        ),
        partial_cens_list=os.path.join(
            config["repeatmasker"]["output_dir"], "status", "all_partial_cens.list"
        ),
    log:
        "logs/repeatmasker/get_complete_correct_cens_bed.log",
    shell:
        """
        grep -f <(cat {input.statuses} | awk '{{
            if ($4 == "false") {{
                split($1, names, ":");
                print names[1]
            }} else {{
                split($1, names, ":");
                print names[1] >> "{output.partial_cens_list}"
            }}
        }}') <(cat {input.beds}) > {output.complete_correct_cens_bed}
        """


use rule extract_and_index_fa as extract_sm_complete_correct_cens with:
    input:
        fa=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}",
            "{sm}_regions.renamed.reort.fa",
        ),
        bed=rules.get_complete_correct_cens_bed.output.complete_correct_cens_bed,
    output:
        seq=os.path.join(
            config["repeatmasker"]["output_dir"],
            "seq",
            "{sm}_complete_correct_cens.fa",
        ),
        idx=os.path.join(
            config["repeatmasker"]["output_dir"],
            "seq",
            "{sm}_complete_correct_cens.fa.fai",
        ),
    params:
        added_cmds="",
    log:
        "logs/repeatmasker/extract_alr_regions_repeatmasker_{sm}.log",


rule merge_all_complete_correct_cens_fa:
    input:
        fa=expand(rules.extract_sm_complete_correct_cens.output.seq, sm=SAMPLE_NAMES),
        idx=expand(rules.extract_sm_complete_correct_cens.output.idx, sm=SAMPLE_NAMES),
    output:
        fa=os.path.join(
            config["repeatmasker"]["output_dir"],
            "seq",
            "all_complete_correct_cens.fa",
        ),
        idx=os.path.join(
            config["repeatmasker"]["output_dir"],
            "seq",
            "all_correct_cens.fa.fai",
        ),
    conda:
        "../env/tools.yaml"
    log:
        "logs/repeatmasker/merge_all_complete_correct_cens.log",
    shell:
        """
        cat {input.fa} > {output.fa}
        samtools faidx {output.fa} 2> {log}
        """


rule remove_partial_cens_rm_out:
    input:
        partial_cens_list=rules.get_complete_correct_cens_bed.output.partial_cens_list,
        rm_out=rules.extract_rm_out_by_chr.output,
    output:
        corrected_rm_out=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "complete_correct_{chr}_cens.fa.out",
        ),
    log:
        "logs/fix_cens_w_repeatmasker/fix_incorrect_mapped_{chr}_cens.log",
    shell:
        """
        grep -v -f {input.partial_cens_list} {input.rm_out} > {output}
        """


rule merge_all_complete_correct_cens_rm:
    input:
        expand(rules.remove_partial_cens_rm_out.output, chr=CHROMOSOMES),
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "all_samples_and_ref_complete_correct_cens.fa.out",
        ),
    shell:
        """
        cat {input} > {output}
        """


rule plot_cens_from_rm_by_chr:
    input:
        script="workflow/scripts/repeatStructure_onlyRM.R",
        rm_out=rules.remove_partial_cens_rm_out.output,
    output:
        repeat_plot_by_chr=os.path.join(
            config["repeatmasker"]["output_dir"],
            "plot",
            "{chr}_cens.pdf",
        ),
    log:
        "logs/fix_cens_w_repeatmasker/plot_{chr}_cens_from_rm.log",
    conda:
        "../env/r.yaml"
    shell:
        """
        Rscript {input.script} {input.rm_out} {output.repeat_plot_by_chr} 2> {log}
        """
