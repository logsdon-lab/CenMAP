include: "common.smk"
include: "utils.smk"


# Check cen status based on repeatmasker annotations of cen and reference.
# Outputs tsv with
# * [original_name, new_name, orientation, is_partial]
rule check_cens_status:
    input:
        rm_out=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "{sm}",
            "{sm}_{chr}_{ctg_name}_correct_ALR_regions.fa.reformatted.out",
        ),
        rm_ref=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "ref",
            "ref_ALR_regions.fa.out",
        ),
    output:
        cens_status=temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "status",
                "{sm}",
                "{sm}_{chr}_{ctg_name}_cens_status.tsv",
            )
        ),
    params:
        edge_len=censtats_status_edge_len,
        edge_perc_alr_thr=censtats_status_edge_perc_alr_thr,
        dst_perc_thr=0.3,
        max_alr_len_thr=censtats_status_max_alr_len_thr,
        # Only allow mapping changes to 13 and 21 if chr13 or chr21.
        restrict_13_21="--restrict_13_21",
    log:
        "logs/fix_cens_w_repeatmasker/check_cens_status_{sm}_{chr}_{ctg_name}.log",
    conda:
        "../envs/py.yaml"
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


def cen_status(wc):
    _ = checkpoints.split_cens_for_rm.get(**wc).output
    fa_glob_pattern = os.path.join(
        config["repeatmasker"]["output_dir"],
        "seq",
        "interm",
        str(wc.sm),
        "{sm}_{chr}_{ctg_name}.fa",
    )
    wcs = glob_wildcards(fa_glob_pattern)

    assert (
        len(wcs.sm) != 0
    ), f"No fasta files found for repeatmasker in {fa_glob_pattern}"
    outputs = []
    for sm, chrom, ctg_name in zip(wcs.sm, wcs.chr, wcs.ctg_name):
        if wc.sm != sm:
            continue
        outputs.extend(
            expand(rules.check_cens_status.output, sm=sm, chr=chrom, ctg_name=ctg_name)
        )
    return outputs


# Generate a key for checking partial and reversed contigs.
# See params.output_format.
rule join_cen_status_and_og_status:
    input:
        status=cen_status,
        rename_key=lambda wc: expand(
            rules.map_collapse_cens.output.renamed_cens_key, sm=wc.sm
        ),
    output:
        all_statuses=temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "status",
                "{sm}_cens_status_no_nucflag.tsv",
            )
        ),
    params:
        # cs - censtats, rn - rename, both
        # (both:new_name_no_ort, cs:is_partial, cs:coords, cs:ort, rn:ort, rn:og_ctg_name)
        output_format="1.1,1.3,1.4,1.2,2.2,2.3",
    log:
        "logs/fix_cens_w_repeatmasker/join_{sm}_cen_status_and_nucflag_status.log",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        join \
        -1 1 -2 1 \
        -a 1 -a 2 \
        -o {params.output_format} \
        <(awk -v OFS="\\t" '{{ split($1, split_name, ":"); print split_name[1],$3,$4,split_name[2]}}' {input.status} | \
            sort -k 1) \
        <(awk -v OFS="\\t" '{{ print $2, ($6 == "false") ? "fwd" : "rev", $1 }}' {input.rename_key} | \
            sort -k 1) > {output} 2> {log}
        """


# Create a final rename key for assembly.
rule get_final_rename_key:
    input:
        statuses=rules.join_cen_status_and_og_status.output,
    output:
        # (original_name, new_name, ort, is_partial, st, end)
        temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "status",
                "{sm}_renamed_cens.tsv",
            )
        ),
    conda:
        "../envs/tools.yaml"
    log:
        "logs/fix_cens_w_repeatmasker/get_{sm}_final_rename_key.log",
    shell:
        """
        {{ awk -v OFS="\\t" '{{
            new_ctg_name=$1;
            og_ctg_name=($6 == "") ? new_ctg_name : $6;
            is_partial=$2;
            split($3, split_coords, "-")
            st=split_coords[1]; end=split_coords[2];
            cs_ort=$4;
            rn_ort=$5;
            if (cs_ort != rn_ort) {{
                print "censtats ort ("cs_ort") and alignment ort ("rn_ort") different for "og_ctg_name >> "{log}"
            }}
            # Rename if reversed.
            if (cs_ort == "rev") {{
                new_ctg_name=gsub("chr", "rc-chr", new_ctg_name)
            }}
            print og_ctg_name,new_ctg_name,cs_ort,is_partial,st,end
        }}' {input.statuses} | \
        sort -k1 ;}} > {output} 2> {log}
        """


# Generate a centromere BED4 file for all cens and mark misassembled cens.
# Also output a repeatmasker reorient key.
rule make_complete_cens_bed:
    input:
        faidx=os.path.join(
            config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fa.fai"
        ),
        final_rename_key=rules.get_final_rename_key.output,
    output:
        # (name, st, end, is_misassembled)
        cen_bed=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            "interm",
            "{sm}_complete_cens.bed",
        ),
        # (old_name, new_name, ort, ctg_len)
        rm_rename_key=os.path.join(
            config["repeatmasker"]["output_dir"],
            "status",
            "{sm}_rm_rename_key.tsv",
        ),
    conda:
        "../envs/tools.yaml"
    log:
        "logs/extract_new_cens_ctgs/fmt_{sm}_new_cens_bed_file.log",
    shell:
        """
        {{ join -1 1 -2 1 <(sort -k1 {input.faidx}) {input.final_rename_key} | \
        awk -v OFS="\\t" '{{
            new_name=$6
            is_partial=$8
            ctg_len=$2
            cen_st=$9
            cen_end=$10
            if ($7 == "rev") {{
                cen_st=ctg_len-cen_end
                cen_end=ctg_len-cen_st
                old_name=gsub("rc-", "", new_name)
                print old_name, new_name, $7, ctg_len >> "{output.rm_rename_key}"
            }}
            print new_name, cen_st, cen_end, is_partial
        }}' ;}} > {output.cen_bed} 2> {log}

        if [ ! -s {output.rm_rename_key} ]; then
            touch {output.rm_rename_key}
        fi
        """


# Generate the original assembly with complete centromeric contigs correctly oriented.
rule rename_reort_asm:
    input:
        fa=os.path.join(config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fa"),
        idx=os.path.join(
            config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fa.fai"
        ),
        rename_key=rules.get_final_rename_key.output,
    output:
        fa=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}_regions.renamed.reort.fa",
        ),
        idx=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}_regions.renamed.reort.fa.fai",
        ),
    params:
        pattern=r"'(\S+)'",
        replacement=lambda wc: "'{kv}'",
    conda:
        "../envs/tools.yaml"
    log:
        "logs/fix_cens_w_repeatmasker/fix_{sm}_asm_orientation.log",
    shell:
        """
        # Get the reverse cens and reverse them.
        # Get all the non-reversed contigs.
        # Then replace the names.
        seqkit replace -p {params.pattern} -r {params.replacement} \
        -k <(cut -f 1,2 {input.rename_key}) \
        <(cat \
            <(seqtk subseq {input.fa} \
                <(awk '$3 == "rev"' {input.rename_key} | cut -f1) | \
                seqtk seq -r) \
            <(seqtk subseq {input.fa} \
                <(grep -v -f <(awk '$3 == "rev"' {input.rename_key} | cut -f1) {input.idx} | cut -f 1)) \
        ) \
        --keep-key > {output.fa} 2> {log}
        samtools faidx {output.fa} 2>> {log}
        """


rule fix_cens_rm_out:
    input:
        rm_rename_key=rules.make_complete_cens_bed.output.rm_rename_key,
        rm_out=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "{sm}_cens.fa.out",
        ),
    output:
        corrected_rm_out=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "reoriented_{sm}_cens.fa.out",
        ),
    log:
        "logs/fix_cens_w_repeatmasker/fix_cens_{sm}_rm_out.log",
    shell:
        """
        # Write everything but the partials and reversed cens.
        grep -v -f \
            <(awk '$3=="rev"' {input.rm_rename_key} | cut -f 1) \
            {input.rm_out} > {output} 2> {log}

        while IFS='' read -r line; do
            original=$(echo "${{line}}" | awk '{{ print $1}}')
            new=$(echo "${{line}}" | awk '{{ print $2}}')
            ctg_len=$(echo "${{line}}" | awk '{{ print $4}}')

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
        done < <(awk -v OFS="\\t" '$3=="rev"' {input.rm_rename_key})
        """


rule merge_fixed_rm_out:
    input:
        expand(rules.fix_cens_rm_out.output, sm=SAMPLE_NAMES),
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "ref",
            "ref_ALR_regions.fa.out",
        ),
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "reoriented_{chr}_cens.fa.out",
        ),
    shell:
        """
        grep -h -P "{wildcards.chr}[:_]" {input} > {output}
        """


use rule plot_rm_out as plot_cens_from_rm_by_chr with:
    input:
        script="workflow/scripts/plot_cens_onlyRM.R",
        rm_out=rules.merge_fixed_rm_out.output,
    output:
        repeat_plot=os.path.join(
            config["repeatmasker"]["output_dir"],
            "plot",
            "{chr}_cens.pdf",
        ),
    log:
        "logs/fix_cens_w_repeatmasker/plot_{chr}_cens_from_rm.log",


rule fix_cens_w_repeatmasker_all:
    input:
        expand(rules.make_complete_cens_bed.output, sm=SAMPLE_NAMES),
        expand(rules.rename_reort_asm.output, sm=SAMPLE_NAMES),
        expand(rules.plot_cens_from_rm_by_chr.output, chr=CHROMOSOMES),
