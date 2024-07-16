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
        edge_len=lambda wc: (
            100_000 if wc.chr in ["chr4", "chr6", "chr14", "chrY"] else 500_000
        ),
        edge_perc_alr_thr=lambda wc: 0.8 if wc.chr in ["chr7"] else 0.95,
        dst_perc_thr=0.3,
        # Edge-case for chrs whose repeats are small and broken up.
        max_alr_len_thr=lambda wc: 0 if wc.chr in ["chrY", "chr11", "chr8"] else 300_000,
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


rule get_cen_corrections_lists:
    input:
        statuses=expand(rules.check_cens_status.output.cens_status, chr=CHROMOSOMES),
    output:
        correct_cens_list=os.path.join(
            config["repeatmasker"]["output_dir"],
            "status",
            "all_complete_correct_cens.list",
        ),
        # (original_name, new_name)
        reverse_cens_key=os.path.join(
            config["repeatmasker"]["output_dir"],
            "status",
            "all_reverse_cens.tsv",
        ),
        partial_cens_list=os.path.join(
            config["repeatmasker"]["output_dir"], "status", "all_partial_cens.list"
        ),
    log:
        "logs/fix_cens_w_repeatmasker/get_complete_correct_cens_bed.log",
    shell:
        """
        cat {input.statuses} | awk -v OFS="\t" '{{
            split($1, names, ":");
            if ($4 == "false") {{
                if ($3 == "fwd") {{
                    print names[1] >> "{output.correct_cens_list}"
                }} else {{
                    # 1. rc-chrX -> rc-rc-chrX (chrX)
                    # 2. rc-rc-chrX -> chrX
                    new_name=names[1]
                    gsub("chr", "rc-chr", new_name)
                    gsub("rc-rc-", "", new_name)
                    print names[1],new_name >> "{output.reverse_cens_key}"
                }}
            }} else {{
                print names[1] >> "{output.partial_cens_list}"
            }}
        }}' 2> {log}
        if [ ! -f {output.partial_cens_list} ]; then
            echo "No partial cens. Creating empty file." > {log}
            touch {output.partial_cens_list}
        fi
        if [ ! -f {output.reverse_cens_key} ]; then
            echo "No cens to reverse. Creating empty file." > {log}
            touch {output.reverse_cens_key}
        fi
        """


rule get_complete_correct_cens_bed:
    input:
        script="workflow/scripts/filter_complete_cens_bed.py",
        bed=rules.get_valid_regions_for_rm.output,
        # Need lengths to determine coordinates for reversed contigs
        fai=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}",
            "{sm}_regions.renamed.reort.fa.fai",
        ),
        correct_cens_list=rules.get_cen_corrections_lists.output.correct_cens_list,
        partial_cens_list=rules.get_cen_corrections_lists.output.partial_cens_list,
        reverse_cens_key=rules.get_cen_corrections_lists.output.reverse_cens_key,
    output:
        complete_correct_cens_bed=os.path.join(
            config["repeatmasker"]["output_dir"],
            "bed",
            "{sm}_complete_correct_ALR_regions.bed",
        ),
    conda:
        "../env/py.yaml"
    log:
        "logs/fix_cens_w_repeatmasker/get_complete_correct_{sm}_cens_bed.log",
    shell:
        """
        python {input.script} \
        -i {input.bed} \
        -l {input.fai} \
        -c {input.correct_cens_list} \
        -p {input.partial_cens_list} \
        -r {input.reverse_cens_key} > {output} 2> {log}
        """


rule fix_ort_asm_final:
    input:
        fa=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}",
            "{sm}_regions.renamed.reort.fa",
        ),
        idx=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}",
            "{sm}_regions.renamed.reort.fa.fai",
        ),
        reverse_cens_key=rules.get_cen_corrections_lists.output.reverse_cens_key,
    output:
        fa=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}",
            "{sm}_regions.renamed.reort.final.fa",
        ),
        idx=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}",
            "{sm}_regions.renamed.reort.final.fa.fai",
        ),
    params:
        tmp_fa="/tmp/{sm}_regions.renamed.reort.final.fa",
    conda:
        "../env/tools.yaml"
    log:
        "logs/fix_cens_w_repeatmasker/fix_{sm }_asm_orientation.log",
    shell:
        """
        # Get the reverse cens and reverse them.
        # Get all the non-reversed contigs.
        seqtk subseq {input.fa} \
            <(grep -f <(cut -f 1 {input.reverse_cens_key}) {input.idx} | cut -f 1) 2> {log} | \
            seqtk seq -r >> {params.tmp_fa}
        seqtk subseq {input.fa} \
            <(grep -v -f <(cut -f 1 {input.reverse_cens_key}) {input.idx} | cut -f 1) 2>> {log} >> {params.tmp_fa}
        # Then replace the names.
        if [ -s {input.reverse_cens_key} ]; then
            seqkit replace -p '(\S+)' -r '{{kv}}' \
            -k {input.reverse_cens_key} {params.tmp_fa} --keep-key > {output.fa} 2> {log}
        else
            mv {params.tmp_fa} {output.fa}
        fi
        samtools faidx {output.fa} 2>> {log}
        """


use rule extract_and_index_fa as extract_sm_complete_correct_cens with:
    input:
        fa=rules.fix_ort_asm_final.output.fa,
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
    log:
        "logs/fix_cens_w_repeatmasker/extract_alr_regions_repeatmasker_{sm}.log",


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
            "all_complete_correct_cens.fa.fai",
        ),
    conda:
        "../env/tools.yaml"
    log:
        "logs/fix_cens_w_repeatmasker/merge_all_complete_correct_cens.log",
    shell:
        """
        cat {input.fa} > {output.fa}
        samtools faidx {output.fa} 2> {log}
        """


rule fix_cens_rm_out:
    input:
        partial_cens_list=rules.get_cen_corrections_lists.output.partial_cens_list,
        reverse_cens_key=rules.get_cen_corrections_lists.output.reverse_cens_key,
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
        # Write everything but the partials and reversed cens.
        grep -v -f \
            <(grep "{wildcards.chr}[_:]" {input.partial_cens_list} | cat - <(grep "{wildcards.chr}[_:]" {input.reverse_cens_key} | cut -f 1)) \
            {input.rm_out} > {output} 2> {log}

        while IFS='' read -r line; do
            original=$(echo "${{line}}" | awk '{{ print $1}}')
            new=$(echo "${{line}}" | awk '{{ print $2}}')
            echo "Replacing ${{original}} with ${{new}} and recalculating coordinates." >> {log}
            # Then replace original name with new name and reverse the output.
            grep "${{original}}" {input.rm_out} | \
                sed "s/${{original}}/${{new}}/g" | \
                tac | \
                awk -v OFS="\\t" '{{
                    # Get start and end coordinates and adjust for reversing.
                    match($5, ":(.+)-", starts);
                    match($5, ".*-(.+)$", ends);
                    new_start=ends[1]-starts[1]-$7+1;
                    new_end=ends[1]-starts[1]-$6+1;
                    $6=new_start;
                    $7=new_end;
                    print \
                }}' >> {output} 2> {log}
        done < <(grep "{wildcards.chr}[_:]" {input.reverse_cens_key})
        """


rule merge_complete_and_correct_rm_out:
    input:
        expand(rules.fix_cens_rm_out.output, chr=CHROMOSOMES),
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
        rm_out=rules.fix_cens_rm_out.output,
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


use rule plot_cens_from_rm_by_chr as plot_cens_from_original_rm_by_chr with:
    input:
        script="workflow/scripts/repeatStructure_onlyRM.R",
        rm_out=rules.extract_rm_out_by_chr.output,
    output:
        repeat_plot_by_chr=os.path.join(
            config["repeatmasker"]["output_dir"],
            "plot",
            "{chr}_cens_original.pdf",
        ),
    log:
        "logs/fix_cens_w_repeatmasker/plot_{chr}_cens_from_original_rm.log",
