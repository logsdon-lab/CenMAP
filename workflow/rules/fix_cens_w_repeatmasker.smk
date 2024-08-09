include: "common.smk"
include: "utils.smk"


rule check_cens_status:
    input:
        rm_out=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "{chr}_cens.fa.out",
        ),
        rm_ref=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "ref",
            "ref_ALR_regions.fa.out",
        ),
    output:
        cens_status=os.path.join(
            config["repeatmasker"]["output_dir"], "status", "{chr}_cens_status.tsv"
        ),
    params:
        edge_len=censtats_status_edge_len,
        edge_perc_alr_thr=censtats_status_edge_perc_alr_thr,
        dst_perc_thr=0.3,
        max_alr_len_thr=censtats_status_max_alr_len_thr,
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


rule join_cen_status_and_nucflag_status:
    input:
        nucflag_statuses=(
            expand(
                os.path.join(config["nucflag"]["output_dir"], "{sm}_cen_status.bed"),
                sm=SAMPLE_NAMES,
            )
            if config.get("nucflag")
            else []
        ),
        statuses=expand(rules.check_cens_status.output.cens_status, chr=CHROMOSOMES),
    output:
        all_statuses=os.path.join(
            config["repeatmasker"]["output_dir"], "status", "cens_status_w_nucflag.tsv"
        ),
    params:
        censtats_base_name_col=5,
        nucflag_base_name_col=1,
        # cs - censtats, nf - nucflag, both
        # [cs:original_name, cs:new_name, cs:ort, cs:is_partial, both:abbreviated_name, nf:status]
        output_format="1.1,1.2,1.3,1.4,0,2.2",
    log:
        "logs/fix_cens_w_repeatmasker/join_cen_status_and_nucflag_status.log",
    shell:
        """
        join \
        -1 {params.censtats_base_name_col} -2 {params.nucflag_base_name_col} \
        -a 1 -a 2 \
        -o {params.output_format} \
        <(awk -v OFS="\\t" '{{ split($1, split_name, ":"); print $0, split_name[1]}}' {input.statuses} | \
            sort -k {params.censtats_base_name_col}) \
        <(cut -f {params.nucflag_base_name_col},4 {input.nucflag_statuses} | \
            sort -k {params.nucflag_base_name_col}) > {output} 2> {log}
        """


rule get_cen_corrections_lists:
    input:
        statuses=(
            rules.join_cen_status_and_nucflag_status.output
            if config.get("nucflag")
            else expand(rules.check_cens_status.output.cens_status, chr=CHROMOSOMES)
        ),
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
        partial_misasm_cens_list=os.path.join(
            config["repeatmasker"]["output_dir"], "status", "all_partial_cens.list"
        ),
    params:
        omit_nucflag="nucflag" not in config,
    log:
        "logs/fix_cens_w_repeatmasker/get_complete_correct_cens_bed.log",
    shell:
        """
        cat {input.statuses} | awk -v OFS="\\t" '{{
            ort=$3; is_partial=$4;
            if ("{params.omit_nucflag}" == "True") {{
                split($1, split_name, ":")
                is_misassembled=""
                name=split_name[1]
            }} else {{
                is_misassembled=$6
                name=$5
            }}
            if (is_partial == "false" && (is_misassembled == "good" || is_misassembled == "")) {{
                if (ort == "fwd") {{
                    print name >> "{output.correct_cens_list}"
                }} else {{
                    # 1. rc-chrX -> rc-rc-chrX (chrX)
                    # 2. rc-rc-chrX -> chrX
                    new_name=name
                    gsub("chr", "rc-chr", new_name)
                    gsub("rc-rc-", "", new_name)
                    print name,new_name >> "{output.reverse_cens_key}"
                }}
            }} else {{
                print name >> "{output.partial_misasm_cens_list}"
            }}
        }}' 2> {log}
        if [ ! -f {output.partial_misasm_cens_list} ]; then
            echo "No partial cens. Creating empty file." > {log}
            touch {output.partial_misasm_cens_list}
        fi
        if [ ! -f {output.reverse_cens_key} ]; then
            echo "No cens to reverse. Creating empty file." > {log}
            touch {output.reverse_cens_key}
        fi
        """


rule get_complete_correct_cens_bed:
    input:
        script="workflow/scripts/filter_complete_cens_bed.py",
        bed=os.path.join(
            config["new_cens"]["output_dir"], "bed", "{sm}_ALR_regions.bed"
        ),
        # Need lengths to determine coordinates for reversed contigs
        fai=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}",
            "{sm}_regions.renamed.reort.fa.fai",
        ),
        correct_cens_list=rules.get_cen_corrections_lists.output.correct_cens_list,
        partial_misasm_cens_list=rules.get_cen_corrections_lists.output.partial_misasm_cens_list,
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
        -p {input.partial_misasm_cens_list} \
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
    conda:
        "../env/tools.yaml"
    log:
        "logs/fix_cens_w_repeatmasker/fix_{sm}_asm_orientation.log",
    shell:
        """
        # Get the reverse cens and reverse them.
        # Get all the non-reversed contigs.
        # Then replace the names.
        if [ -s {input.reverse_cens_key} ]; then
            seqkit replace -p '(\S+)' -r '{{kv}}' \
            -k {input.reverse_cens_key} \
            <(cat \
                <(seqtk subseq {input.fa} \
                    <(grep -f <(cut -f 1 {input.reverse_cens_key}) {input.idx} | cut -f 1) | \
                    seqtk seq -r) \
                <(seqtk subseq {input.fa} \
                    <(grep -v -f <(cut -f 1 {input.reverse_cens_key}) {input.idx} | cut -f 1)) \
            ) \
            --keep-key > {output.fa} 2> {log}
        else
            cp {input.fa} {output.fa}
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


rule get_reverse_cens_key_with_ctg_len:
    input:
        reverse_cens_key=rules.get_cen_corrections_lists.output.reverse_cens_key,
        fais=expand(
            os.path.join(
                config["concat_asm"]["output_dir"],
                "{sm}",
                "{sm}_regions.renamed.reort.fa.fai",
            ),
            sm=SAMPLE_NAMES,
        ),
    output:
        reverse_cens_key_w_len=os.path.join(
            config["repeatmasker"]["output_dir"],
            "status",
            "all_reverse_cens_w_ctg_len.tsv",
        ),
    log:
        "logs/get_reverse_cens_key_with_ctg_len.log",
    shell:
        """
        {{ join <(sort -k1 {input.reverse_cens_key}) <(sort -k1 {input.fais}) | \
        sed 's/ /\\t/g'| \
        cut -f 1,2,3 ;}} > {output} 2> {log}
        """


rule fix_cens_rm_out:
    input:
        partial_misasm_cens_list=rules.get_cen_corrections_lists.output.partial_misasm_cens_list,
        reverse_cens_key_w_len=rules.get_reverse_cens_key_with_ctg_len.output,
        rm_out=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "{chr}_cens.fa.out",
        ),
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
            <(grep "{wildcards.chr}[_:]" {input.partial_misasm_cens_list} | cat - <(grep "{wildcards.chr}[_:]" {input.reverse_cens_key_w_len} | cut -f 1)) \
            {input.rm_out} > {output} 2> {log}

        while IFS='' read -r line; do
            original=$(echo "${{line}}" | awk '{{ print $1}}')
            new=$(echo "${{line}}" | awk '{{ print $2}}')
            ctg_len=$(echo "${{line}}" | awk '{{ print $3}}')

            echo "Replacing ${{original}} with ${{new}} and recalculating coordinates." >> {log}
            # Then replace original name with new name and reverse the output.
            # Also include contig len to account for reversal.
            grep "${{original}}" {input.rm_out} | \
                sed "s/${{original}}/${{new}}/g" | \
                tac | \
                awk -v CTG_LEN="${{ctg_len}}" -v OFS="\\t" '{{
                    # Get start and end coordinates and adjust for reversing.
                    match($5, "(.+):", ctgs);
                    match($5, ":(.+)-", starts);
                    match($5, ".*-(.+)$", ends);
                    new_start=ends[1]-starts[1]-$7+1;
                    new_end=ends[1]-starts[1]-$6+1;
                    $6=new_start;
                    $7=new_end;
                    # Then rename ctg.
                    new_ctg_start=CTG_LEN-ends[1]+1;
                    new_ctg_end=CTG_LEN-starts[1]+1;
                    new_ctg_name=ctgs[1]":"new_ctg_start"-"new_ctg_end;
                    $5=new_ctg_name;
                    print \
                }}' >> {output} 2> {log}
        done < <(grep "{wildcards.chr}[_:]" {input.reverse_cens_key_w_len})
        """


use rule plot_rm_out as plot_cens_from_rm_by_chr with:
    input:
        script="workflow/scripts/repeatStructure_onlyRM.R",
        rm_out=rules.fix_cens_rm_out.output,
    output:
        repeat_plot=os.path.join(
            config["repeatmasker"]["output_dir"],
            "plot",
            "{chr}_cens.pdf",
        ),
    log:
        "logs/fix_cens_w_repeatmasker/plot_{chr}_cens_from_rm.log",


rule fix_cens_w_repeatmasker_only:
    input:
        expand(rules.fix_ort_asm_final.output, sm=SAMPLE_NAMES),
        expand(rules.extract_sm_complete_correct_cens.output, sm=SAMPLE_NAMES),
        expand(rules.merge_all_complete_correct_cens_fa.output, sm=SAMPLE_NAMES),
        expand(rules.plot_cens_from_rm_by_chr.output, chr=CHROMOSOMES),
