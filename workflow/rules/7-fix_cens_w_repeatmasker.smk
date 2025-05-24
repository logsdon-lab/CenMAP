include: "common.smk"
include: "utils.smk"
include: "1-concat_asm.smk"
include: "4-ident_cen_ctgs.smk"
include: "6-repeatmasker.smk"


FIX_RM_OUTDIR = join(OUTPUT_DIR, "7-fix_cens_w_repeatmasker")
FIX_RM_LOGDIR = join(LOG_DIR, "7-fix_cens_w_repeatmasker")
FIX_RM_BMKDIR = join(BMK_DIR, "7-fix_cens_w_repeatmasker")


# Check cen status based on repeatmasker annotations of cen and reference.
# Outputs tsv with
# * [original_name, new_name, orientation, is_partial]
rule check_cens_status:
    input:
        rm_out=rules.reformat_repeatmasker_output.output,
        rm_ref=rules.merge_control_repeatmasker_output.output,
    output:
        cens_status=join(
            FIX_RM_OUTDIR,
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
        join(FIX_RM_LOGDIR, "check_cens_status_{sm}_{fname}.log"),
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
    outdir = checkpoints.split_cens_for_rm.get(**wc).output[0]
    fa_glob_pattern = join(outdir, "{fname}.fa")
    wcs = glob_wildcards(fa_glob_pattern)
    fnames = wcs.fname
    return expand(rules.check_cens_status.output, sm=wc.sm, fname=fnames)


# Generate a key for checking partial and reversed contigs.
rule make_complete_cens_bed:
    input:
        script=workflow.source_path("../scripts/make_complete_cens_bed.py"),
        status=cen_status,
        faidx=rules.concat_asm.output.idx,
    output:
        # (new_name, st, end, is_misassembled, ctg_name, ctg_len)
        cen_bed=join(
            FIX_RM_OUTDIR,
            "bed",
            "interm",
            "{sm}_complete_cens.bed",
        ),
    params:
        use_censtats_new_name=(
            "--use_new_name"
            if config["repeatmasker"].get("use_censtats_new_name")
            else ""
        )
    log:
        join(FIX_RM_LOGDIR, "get_complete_cens_bed_{sm}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} -i {input.status} -f {input.faidx} {params.use_censtats_new_name} > {output} 2> {log}
        """


# Generate the original assembly with complete centromeric contigs correctly oriented.
rule rename_reort_asm:
    input:
        fa=rules.concat_asm.output.fa,
        idx=rules.concat_asm.output.idx,
        cens_bed=ancient(rules.make_complete_cens_bed.output),
    output:
        fa=join(
            CONCAT_ASM_OUTDIR,
            "{sm}-asm-renamed-reort.fa",
        ),
        idx=join(
            CONCAT_ASM_OUTDIR,
            "{sm}-asm-renamed-reort.fa.fai",
        ),
    params:
        pattern=r"'^(\S+)\s*'",
        replacement=lambda wc: "'{kv}'",
    conda:
        "../envs/tools.yaml"
    log:
        join(FIX_RM_LOGDIR, "fix_{sm}_asm_orientation.log"),
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


# Rename and reorient original RM annotations with changes.
rule fix_cens_rm_out:
    input:
        script=workflow.source_path("../scripts/rename_reorient_rm_out.py"),
        rm_rename_key=ancient(rules.make_complete_cens_bed.output),
        rm_out=ancient(rules.format_repeatmasker_output.output),
    output:
        corrected_rm_out=join(
            FIX_RM_OUTDIR,
            "repeats",
            "all",
            "reoriented_{sm}_cens.fa.out",
        ),
    log:
        join(FIX_RM_LOGDIR, "fix_cens_{sm}_rm_out.log"),
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} -i {input.rm_out} -k {input.rm_rename_key} > {output} 2> {log}
        """


rule create_fixed_rm_bed:
    input:
        script=workflow.source_path("../scripts/create_rm_bed.py"),
        rm_out=[
            expand(rules.fix_cens_rm_out.output, sm=SAMPLE_NAMES),
            rules.merge_control_repeatmasker_output.output,
        ],
    output:
        rm_bed=join(
            FIX_RM_OUTDIR,
            "bed",
            "{chr}",
            "rm.bed",
        ),
    params:
        chr_rgx="{chr}[:_]",
        color_mapping=config["repeatmasker"]["repeat_colors"],
    log:
        join(FIX_RM_LOGDIR, "create_fixed_rm_bed_{chr}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        shell_create_rm_bed


rule modify_rm_cenplot_tracks:
    input:
        plot_layout=workflow.source_path("../scripts/cenplot_repeatmasker_plot.toml"),
        infiles=rules.create_fixed_rm_bed.output,
    output:
        plot_layout=join(
            FIX_RM_OUTDIR,
            "plots",
            "{chr}_cens_fixed.yaml",
        ),
    run:
        format_toml_path(
            input_plot_layout=input.plot_layout,
            output_plot_layout=output.plot_layout,
            indir=os.path.abspath(os.path.dirname(str(input.infiles[0]))),
            **dict(params.items()),
        )


rule plot_fixed_rm_bed_by_chr:
    input:
        bed_files=[rules.create_fixed_rm_bed.output],
        script=workflow.source_path("../scripts/plot_multiple_cen.py"),
        plot_layout=rules.modify_rm_cenplot_tracks.output,
    output:
        plots=multiext(
            join(
                FIX_RM_OUTDIR,
                "plots",
                "{chr}_cens",
            ),
            ".pdf",
            ".png",
        ),
        plot_dir=directory(
            join(
                FIX_RM_OUTDIR,
                "plots",
                "{chr}_cens",
            )
        ),
    threads: 4
    log:
        join(FIX_RM_LOGDIR, "plot_fixed_rm_bed_{chr}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        shell_plot_multiple_cen


rule fix_cens_w_repeatmasker_all:
    input:
        expand(rules.rename_reort_asm.output, sm=SAMPLE_NAMES),
        expand(rules.plot_fixed_rm_bed_by_chr.output, chr=CHROMOSOMES),
    default_target: True
