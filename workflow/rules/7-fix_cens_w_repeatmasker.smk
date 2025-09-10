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
        entropy_bed=rules.calculate_entropy.output,
        rm_out=rules.reformat_repeatmasker_output.output,
    output:
        # (chrom, st, end, old_chrom)
        bed=temp(
            join(
                FIX_RM_OUTDIR,
                "entropy",
                "interm",
                "{sm}_{fname}.bed",
            )
        ),
    log:
        join(FIX_RM_LOGDIR, "filter_entropy_bed_{sm}_{fname}.log"),
    params:
        script=workflow.source_path("../scripts/filter_entropy_bed.py"),
        trim_to_repeats=(
            f"--trim_to_repeats {' '.join(config['repeatmasker']['trim_to_repeats'])}"
            if config["repeatmasker"]["trim_to_repeats"]
            else ""
        ),
        bp_merge=config["ident_cen_ctgs"]["bp_merge"],
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {params.script} \
        -i {input.entropy_bed} \
        -r {input.rm_out} \
        {params.trim_to_repeats} \
        -d {params.bp_merge} > {output} 2> {log}
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
        # (sm_ctg, st, end, old_sm_ctg)
        beds=ancient(valid_beds_by_cen_entropy),
        # (ctg, sm_ctg, ctg_len)
        rename_key=rules.map_chroms.output,
        idx=rules.rename_reort_asm.output.idx,
    output:
        # (sm_ctg, st, end, ctg, ctg_len)
        cen_bed=join(
            FIX_RM_OUTDIR,
            "bed",
            "interm",
            "{sm}_complete_cens.bed",
        ),
        # (old_sm_ctg, sm_ctg)
        rename_key=temp(
            join(
                FIX_RM_OUTDIR,
                "bed",
                "interm",
                "{sm}_complete_cens_rename_key.tsv",
            )
        ),
    conda:
        "../envs/tools.yaml"
    params:
        bp_slop=config["ident_cen_ctgs"]["bp_slop"],
    shell:
        """
        sort -k1,1 -k 2,2n {input.beds}  | \
        join -1 1 -2 2 - <(sort -k2,2 {input.rename_key}) | \
        awk -v OFS="\\t" '{{ $1=$1; print }}' | \
        bedtools slop -i - -g {input.idx} -b {params.bp_slop} | \
        awk -v OFS="\\t" '{{
            # Adjust for seqtk 1-index
            $2=($2 == 0) ? 1 : $2
            new_name=$1":"$2+1"-"$3
            print $4, new_name >> "{output.rename_key}"
            print $1, $2, $3, $5, $6
        }}' > {output.cen_bed}
        """


# Filter original RM annotations.
rule fix_cens_rm_out:
    input:
        bed=rules.make_complete_cens_bed.output.cen_bed,
        rename_key=rules.make_complete_cens_bed.output.rename_key,
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
        {{ awk -v OFS="\\t" '{{
            if (FNR == NR) {{ kv[$1]=$2; next; }};
            new_name=kv[$5];
            match($5, ":(.+)-", ctg_st);
            match($5, "^(.+):", ctg_name);
            # Convert to absolute coordinates
            $6=$6+ctg_st[1];
            $7=$7+ctg_st[1];
            # Make bed-like
            if (new_name) {{
                $5=new_name;
                print ctg_name[1], $6, $7, $0
            }}
        }}' {input.rename_key} <(cut -f1-15 {input.rm_out}) | \
        bedtools intersect -a - -b {input.bed} | \
        cut -f 4-18 ;}} > {output} 2> {log}
        """


rule create_fixed_rm_bed:
    input:
        rm_out=expand(rules.fix_cens_rm_out.output, sm=SAMPLE_NAMES),
    output:
        rm_bed=join(
            FIX_RM_OUTDIR,
            "bed",
            "{chr}",
            "rm.bed",
        ),
    params:
        chr_rgx="{chr}[:_-]",
        color_mapping=config["repeatmasker"]["repeat_colors"],
        script=workflow.source_path("../scripts/create_rm_bed.py"),
        to_abs="",
    log:
        join(FIX_RM_LOGDIR, "create_fixed_rm_bed_{chr}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        shell_create_rm_bed


rule plot_fixed_rm_bed_by_chr:
    input:
        rm=[
            (
                rules.create_ref_rm_bed.output
                if config["repeatmasker"]["ref_repeatmasker_output"]
                else []
            ),
            rules.create_fixed_rm_bed.output,
        ],
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
    params:
        script=workflow.source_path("../scripts/plot_multiple_cen.py"),
        plot_layout=workflow.source_path("../scripts/cenplot_repeatmasker_plot.yaml"),
        output_prefix=lambda wc, output: os.path.splitext(output.plots[0])[0],
        json_file_str=lambda wc, input: json.dumps(dict(input)),
        options="",
        omit_if_empty="",
        ref_ax_idx="--ref_ax_idx 0",
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
