include: "common.smk"
include: "utils.smk"
include: "5-ident_cen_ctgs.smk"
include: "7-fix_cens_w_repeatmasker.smk"


# Get sample names with subdirs in cdr_finder.input_bam_dir
SAMPLE_NAMES_BAM = set(
    glob_wildcards(join(config["cdr_finder"]["input_bam_dir"], "{sm}")).sm
)
SAMPLE_NAMES_INTERSECTION = SAMPLE_NAMES_BAM.intersection(set(SAMPLE_NAMES))
print(
    f"Using {len(SAMPLE_NAMES_INTERSECTION)} out of {len(SAMPLE_NAMES)} samples "
    f"that have subdirs in {config['cdr_finder']['input_bam_dir']} for CDR-Finder.",
    file=sys.stderr,
)
ALIGNER = config["cdr_finder"].get("aligner", "minimap2")
DEF_ALIGNER_OPTS = {
    "minimap2": "-y -a --eqx --cs -x lr:hqae -I8g -s 4000",
    # https://www.biorxiv.org/content/10.1101/2024.11.01.621587v1
    "winnowmap": "-y -a --eqx --cs -x map-ont -I8g -s 4000",
}
try:
    ALIGNER_OPTS = config["cdr_finder"].get("aligner_opts", DEF_ALIGNER_OPTS[ALIGNER])
except KeyError:
    raise ValueError(
        f"Invalid aligner option ({ALIGNER}) for CDR-Finder. Choose: {tuple(DEF_ALIGNER_OPTS.keys())} "
    )


CDR_FINDER_OUTDIR = join(OUTPUT_DIR, "8-cdr_finder")
CDR_FINDER_LOGDIR = join(LOG_DIR, "8-cdr_finder")
CDR_FINDER_BMKDIR = join(BMK_DIR, "8-cdr_finder")


wildcard_constraints:
    sm="|".join(SAMPLE_NAMES_INTERSECTION),


CDR_ALIGN_CFG = {
    "samples": [
        {
            "name": sm,
            "asm_fa": expand(rules.create_final_asm.output.fa, sm=sm)[0],
            "read_dir": join(config["cdr_finder"]["input_bam_dir"], sm),
            "read_rgx": config["cdr_finder"]["bam_rgx"],
        }
        for sm in SAMPLE_NAMES_INTERSECTION
    ],
    "aligner": ALIGNER,
    "aligner_opts": ALIGNER_OPTS,
    "mem_aln": config["cdr_finder"]["aln_mem"],
    "threads_aln": config["cdr_finder"]["aln_threads"],
    "benchmarks_dir": CDR_FINDER_BMKDIR,
    "logs_dir": CDR_FINDER_LOGDIR,
    "output_dir": join(CDR_FINDER_OUTDIR, "aln"),
    "samtools_view_flag": 2308,
}


module CDR_Align:
    snakefile:
        "Snakemake-Aligner/workflow/Snakefile"
    config:
        CDR_ALIGN_CFG


use rule * from CDR_Align as cdr_aln_*


# Avoid running RM again.
rule cdr_make_repeatmasker_bed:
    input:
        rules.fix_cens_rm_out.output,
    output:
        # [chrom_w_no_coords, st_abs, end_abs, repeat_type]
        rm_bed=temp(
            join(
                CDR_FINDER_OUTDIR,
                "bed",
                "{sm}_rm.bed",
            )
        ),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        awk -v OFS="\\t" '{{
            match($5, "^(.+):", chrom);
            print chrom[1], $6, $7, $10
        }}' {input} > {output}
        """


# Use SRF regions as putative alpha-satellite regions.
rule cdr_make_srf_putative_alr_regions:
    input:
        bed=rules.extract_filter_monomers.output.bed_mon,
        rename_key=rules.create_rename_key.output if CHROMOSOMES else [],
    output:
        join(
            CDR_FINDER_OUTDIR,
            "bed",
            "{sm}_region.bed",
        ),
    params:
        alr_threshold=config.get("cdr_finder", {}).get("alr_threshold"),
        bp_alr_merge=config.get("cdr_finder", {}).get("bp_alr_merge"),
        format_cmd=lambda wc, input: (
            f"""join - <(sort -k1,1 {input.rename_key}) | \             awk -v OFS="\\t" '{{{{                 if ($4 ~ "rc-chr") {{{{ st=$5-$3; end=$5-$2; }}}} else {{{{ st=$2; end=$3; }}}};                 print $4, st, end, "ALR/Alpha"             }}}}' | sort | uniq"""
            if input.rename_key
            else f"""awk -v OFS="\\t" '{{{{ $1="{wc.sm}_"$1; print $1, $2, $3, "ALR/Alpha" }}}}'"""
        ),
    log:
        join(CDR_FINDER_LOGDIR, "create_srf_cdr_finder_regions_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ srf-n-trf regions \
        -b {input.bed} \
        -d {params.bp_alr_merge} \
        -m {params.alr_threshold} | \
        sort -k1,1 | \
        {params.format_cmd} ;}} > {output} 2> {log}
        """


# Pass CDR config here.
CDR_FINDER_CONFIG = {
    "output_dir": CDR_FINDER_OUTDIR,
    "log_dir": CDR_FINDER_LOGDIR,
    "benchmark_dir": CDR_FINDER_BMKDIR,
    "restrict_alr": True,
    "samples": {
        sm: {
            "fasta": expand(rules.create_final_asm.output.fa, sm=sm),
            "regions": expand(rules.make_complete_cens_bed.output.cen_bed, sm=sm),
            "bam": expand(
                rules.cdr_aln_merge_read_asm_alignments.output.alignment, sm=sm
            ),
            "repeatmasker": expand(
                (
                    rules.cdr_make_srf_putative_alr_regions.output
                    if config["cdr_finder"].get("use_srf_regions")
                    else rules.cdr_make_repeatmasker_bed.output.rm_bed
                ),
                sm=sm,
            ),
        }
        for sm in SAMPLE_NAMES_INTERSECTION
    },
    **config["cdr_finder"],
}


module CDR_Finder:
    snakefile:
        "CDR-Finder/workflow/Snakefile"
    config:
        CDR_FINDER_CONFIG


use rule * from CDR_Finder as cdr_*


rule merge_cdr_beds:
    input:
        bed=expand(
            rules.make_complete_cens_bed.output.cen_bed, sm=SAMPLE_NAMES_INTERSECTION
        ),
        cdr_output=expand(rules.cdr_call_cdrs.output, sample=SAMPLE_NAMES_INTERSECTION),
        methyl_cdr_output=expand(
            rules.cdr_aggregate_bed_files.output.methyl_bed_all,
            sample=SAMPLE_NAMES_INTERSECTION,
        ),
    output:
        reorient_cdr_output=join(
            CDR_FINDER_OUTDIR,
            "bed",
            "all_cdrs.bed",
        ),
        reorient_methyl_cdr_output=join(
            CDR_FINDER_OUTDIR,
            "bed",
            "all_binned_freq.bed",
        ),
    shell:
        """
        cut -f1-3 {input.cdr_output} | \
        sort -k1,1 | \
        join - <(cat {input.bed} | sort -k1,1) | \
        awk -v OFS="\\t" '{{
            if ($2 >= $4 && $3 <= $5) {{
                print $1":"$4+1"-"$5, $2,$3
            }}
        }}' > {output.reorient_cdr_output}
        sort -k1,1 {input.methyl_cdr_output}  | \
        join - <(cat {input.bed} | sort -k1,1) | \
        awk -v OFS="\\t" '{{
            if ($2 >= $6 && $3 <= $7) {{
                print $1":"$6+1"-"$7,$2,$3,$4,$5
            }}
        }}' > {output.reorient_methyl_cdr_output}
        """


rule cdr_finder_all:
    input:
        expand(rules.cdr_cdr_plots.output, sample=SAMPLE_NAMES_INTERSECTION),
        rules.merge_cdr_beds.output,
    default_target: True
