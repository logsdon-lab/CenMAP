include: "common.smk"
include: "utils.smk"
include: "5-ident_cen_ctgs.smk"
include: "6-repeatmasker.smk"


FIX_RM_OUTDIR = join(OUTPUT_DIR, "7-fix_cens_w_repeatmasker")
FIX_RM_LOGDIR = join(LOG_DIR, "7-fix_cens_w_repeatmasker")
FIX_RM_BMKDIR = join(BMK_DIR, "7-fix_cens_w_repeatmasker")


rule calculate_entropy:
    input:
        rm_out=rules.reformat_repeatmasker_output.output,
    output:
        # BED9
        entropy_bed=join(
            FIX_RM_OUTDIR,
            "entropy",
            "{sm}",
            "{fname}.bed",
        ),
    params:
        outdir=lambda wc, output: os.path.dirname(output[0]),
        window=config["repeatmasker"]["bp_shannon_window"],
    log:
        join(FIX_RM_LOGDIR, "calculate_entropy_{sm}_{fname}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        """
        censtats entropy \
        -i <(awk -v OFS="\\t" '{{ print $5, $6, $7, $10, $9 }}' {input.rm_out}) \
        -w {params.window} \
        -o {params.outdir} 2> {log}
        """


# Filter valid cens based on entropy and overlap with ALR/Alpha
# Then trim.
rule filter_entropy_bed:
    input:
        script=workflow.source_path("../scripts/filter_entropy_bed.py"),
        entropy_bed=rules.calculate_entropy.output,
        rm_out=rules.reformat_repeatmasker_output.output,
    output:
        # (chrom, st, end)
        bed=join(
            FIX_RM_OUTDIR,
            "bed",
            "{sm}",
            "{fname}.bed",
        ),
    log:
        join(FIX_RM_LOGDIR, "filter_entropy_bed_{sm}_{fname}.log"),
    params:
        trim_to_repeats=(
            f"--trim_to_repeats {' '.join(config['repeatmasker']['trim_to_repeats'])}"
            if config["repeatmasker"]["trim_to_repeats"]
            else ""
        ),
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} \
        -i {input.entropy_bed} \
        -r {input.rm_out} \
        {params.trim_to_repeats} > {output} 2> {log}
        """


def valid_beds_by_cen_entropy(wc):
    outdir = checkpoints.split_cens_for_rm.get(**wc).output[0]
    fa_glob_pattern = join(outdir, "{fname}.fa")
    wcs = glob_wildcards(fa_glob_pattern)
    fnames = wcs.fname
    return expand(rules.filter_entropy_bed.output, sm=wc.sm, fname=fnames)


# Merge valid cens coordinates.
rule make_complete_cens_bed:
    input:
        beds=valid_beds_by_cen_entropy,
        rename_key=rules.map_chroms.output,
        idx=rules.rename_reort_asm.output.idx,
    output:
        # (new_name, st, end, ctg, ctg_len)
        cen_bed=join(
            FIX_RM_OUTDIR,
            "bed",
            "interm",
            "{sm}_complete_cens.bed",
        ),
    conda:
        "../envs/tools.yaml"
    params:
        bp_slop=config["ident_cen_ctgs"]["bp_slop"],
    shell:
        """
        cat {input.beds} | \
        sort -k1,1 -k 2,2n | \
        join -1 1 -2 2 - <(sort -k2,2 {input.rename_key}) | \
        awk -v OFS="\\t" '{{$1=$1; print}}' | \
        bedtools slop -i - -g {input.idx} -b {params.bp_slop} > {output}
        """


# Filter original RM annotations.
rule fix_cens_rm_out:
    input:
        bed=rules.make_complete_cens_bed.output,
        rm_out=rules.format_repeatmasker_output.output,
    output:
        corrected_rm_out=join(
            FIX_RM_OUTDIR,
            "repeats",
            "all",
            "{sm}_cens.fa.out",
        ),
    log:
        join(FIX_RM_LOGDIR, "fix_cens_{sm}_rm_out.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        grep -f <(awk '{{print $1":"$2"-"$3}}' {input.bed}) {input.rm_out} > {output} 2> {log}
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
        script=workflow.source_path("../scripts/format_cenplot_toml.py"),
        plot_layout=workflow.source_path("../scripts/cenplot_repeatmasker_plot.toml"),
        infiles=rules.create_fixed_rm_bed.output,
    output:
        plot_layout=join(
            FIX_RM_OUTDIR,
            "plots",
            "{chr}_cens_fixed.yaml",
        ),
    conda:
        "../envs/py.yaml"
    log:
        join(FIX_RM_LOGDIR, "modify_rm_cenplot_tracks_{chr}.log"),
    params:
        indir=lambda wc, input: os.path.abspath(os.path.dirname(str(input.infiles[0]))),
    shell:
        """
        python {input.script} \
        -i {input.plot_layout} \
        -o {output.plot_layout} \
        -k indir={params.indir} &> {log}
        """


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
        expand(rules.plot_fixed_rm_bed_by_chr.output, chr=CHROMOSOMES),
    default_target: True
