include: "common.smk"
include: "utils.smk"


# Check cen status based on repeatmasker annotations of cen and reference.
# Outputs tsv with
# * [original_name, new_name, orientation, is_partial]
rule check_cens_status:
    input:
        rm_out=os.path.join(
            config["repeatmasker"]["output_dir"], "repeats", "{sm}", "{fname}.fa.out"
        ),
        rm_ref=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "ref",
            "ref_ALR_regions.fa.out",
        ),
    output:
        cens_status=os.path.join(
            config["repeatmasker"]["output_dir"],
            "status",
            "{sm}",
            "{fname}_cens_status.tsv",
        ),
    params:
        edge_len=lambda wc: censtats_status_edge_len(get_chrom_name(wc.fname)),
        edge_perc_alr_thr=lambda wc: censtats_status_edge_perc_alr_thr(
            get_chrom_name(wc.fname)
        ),
        dst_perc_thr=0.3,
        max_alr_len_thr=lambda wc: censtats_status_max_alr_len_thr(
            get_chrom_name(wc.fname)
        ),
        # Only allow mapping changes to 13 and 21 if chr13 or chr21.
        restrict_13_21="--restrict_13_21",
    log:
        "logs/fix_cens_w_repeatmasker/check_cens_status_{sm}_{fname}.log",
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
    try:
        _ = checkpoints.split_cens_for_rm.get(**wc).output
    except AttributeError:
        pass
    wcs = glob_wildcards(
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "seq",
            "interm",
            str(wc.sm),
            "{fname}.fa",
        )
    )
    fnames = wcs.fname
    assert (
        len(fnames) != 0
    ), f"No fasta files found for repeatmasker in {fa_glob_pattern}"

    return expand(rules.check_cens_status.output, sm=wc.sm, fname=fnames)


# Generate a key for checking partial and reversed contigs.
rule make_complete_cens_bed:
    input:
        script="workflow/scripts/make_complete_cens_bed.py",
        status=cen_status,
        faidx=os.path.join(
            config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fa.fai"
        ),
    output:
        # (new_name, st, end, is_misassembled, ctg_name, ctg_len)
        cen_bed=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            "interm",
            "{sm}_complete_cens.bed",
        ),
    log:
        "logs/fix_cens_w_repeatmasker/get_complete_cens_bed_{sm}.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} -i {input.status} -f {input.faidx} > {output} 2> {log}
        """


# Generate the original assembly with complete centromeric contigs correctly oriented.
rule rename_reort_asm:
    input:
        fa=os.path.join(config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fa"),
        idx=os.path.join(
            config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fa.fai"
        ),
        cens_bed=rules.make_complete_cens_bed.output,
    output:
        fa=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}-asm-renamed-reort.fa",
        ),
        idx=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}-asm-renamed-reort.fa.fai",
        ),
    params:
        pattern=r"'^(\S+)\s*'",
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
        -k <(awk -v OFS="\\t" '{{ print $5, $1}}' {input.cens_bed}) \
        <(cat \
            <(seqtk subseq {input.fa} \
                <(awk '$1 ~ "rc-chr"' {input.cens_bed} | cut -f5) | \
                seqtk seq -r) \
            <(seqtk subseq {input.fa} \
                <(grep -v -f <(awk '$1 ~ "rc-chr"' {input.cens_bed} | cut -f5) {input.idx} | cut -f 1)) \
        ) \
        --keep-key > {output.fa} 2> {log}
        samtools faidx {output.fa} 2>> {log}
        """


rule fix_cens_rm_out:
    input:
        script="workflow/scripts/rename_reorient_rm_out.py",
        rm_rename_key=rules.make_complete_cens_bed.output,
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
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} -i {input.rm_out} -k {input.rm_rename_key} > {output} 2> {log}
        """


use rule create_rm_bed as create_fixed_rm_bed with:
    input:
        script="workflow/scripts/create_rm_bed.py",
        # Absolute coordinates
        rm_out=[
            expand(rules.fix_cens_rm_out.output, sm=SAMPLE_NAMES),
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "repeats",
                "ref",
                "ref_ALR_regions.fa.abs.out",
            ),
        ],
    output:
        rm_bed=os.path.join(
            config["repeatmasker"]["output_dir"],
            "bed",
            "{chr}",
            "rm.bed",
        ),
    params:
        chr_rgx="{chr}[:_]",
        color_mapping=config["repeatmasker"]["repeat_colors"],
    log:
        "logs/fix_cens_w_repeatmasker/create_fixed_rm_bed_{chr}.log",


use rule plot_multiple_cen as plot_fixed_rm_bed_by_chr with:
    input:
        bed_files=[rules.create_fixed_rm_bed.output],
        script="workflow/scripts/plot_multiple_cen.py",
        # Fixed tracks. No differnce between this and og.
        plot_layout=expand(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "plots",
                "{chr}_cens_{typ}.yaml",
            ),
            chr="{chr}",
            typ="fixed",
        ),
    output:
        plots=multiext(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "plots",
                "{chr}_cens",
            ),
            ".pdf",
            ".png",
        ),
        plot_dir=directory(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "plots",
                "{chr}_cens",
            )
        ),
    log:
        "logs/fix_cens_w_repeatmasker/plot_fixed_rm_bed_{chr}.log",


rule fix_cens_w_repeatmasker_all:
    input:
        expand(rules.make_complete_cens_bed.output, sm=SAMPLE_NAMES),
        expand(rules.rename_reort_asm.output, sm=SAMPLE_NAMES),
        expand(rules.plot_fixed_rm_bed_by_chr.output, chr=CHROMOSOMES),
    default_target: True
