include: "common.smk"
include: "utils.smk"
include: "1-concat_asm.smk"
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
            "asm_fa": expand(rules.concat_asm.output.fa, sm=sm),
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
        github(
            "logsdon-lab/Snakemake-Aligner", path="workflow/Snakefile", branch="main"
        )
    config:
        CDR_ALIGN_CFG


use rule * from CDR_Align as cdr_aln_*


rule get_original_coords:
    input:
        # (sample_chr_ctg, st, end, is-misassembled)
        bed=ancient(rules.make_complete_cens_bed.output),
        asm_faidx=rules.concat_asm.output.idx,
    output:
        og_coords=join(
            CDR_FINDER_OUTDIR,
            "bed",
            "{sm}_complete_correct_ALR_regions.og_coords.bed",
        ),
        og_coords_key=join(
            CDR_FINDER_OUTDIR,
            "bed",
            "{sm}_complete_correct_ALR_regions.og_coords.key.bed",
        ),
    log:
        join(CDR_FINDER_LOGDIR, "get_original_coords_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ join -1 1 -2 1 \
            <(sort -k 5 {input.bed} | awk -v OFS="\\t" '{{print $5, $2, $3, $1}}') \
            <(cut -f 1,2 {input.asm_faidx} | sort -k 1) | \
        awk -v OFS="\\t" '{{
            ctg=$1
            final_ctg=$4
            len=$5
            is_reversed=($4 ~ "rc")
            start=$2; end=$3;
            adj_start=$2; adj_end=$3;
            if (is_reversed) {{
                adj_start=len-end
                adj_end=len-start
            }}
            print ctg, adj_start, adj_end > "{output.og_coords}"
            print ctg, adj_start, adj_end, final_ctg, start, end, len
        }}';}} > {output.og_coords_key} 2> {log}
        """


# Pass CDR config here.
CDR_FINDER_CONFIG = {
    "output_dir": CDR_FINDER_OUTDIR,
    "log_dir": CDR_FINDER_LOGDIR,
    "benchmark_dir": CDR_FINDER_BMKDIR,
    "restrict_alr": True,
    "samples": {
        sm: {
            "fasta": expand(rules.concat_asm.output.fa, sm=sm),
            "regions": expand(rules.get_original_coords.output.og_coords, sm=sm),
            "bam": expand(
                rules.cdr_aln_merge_read_asm_alignments.output.alignment, sm=sm
            ),
        }
        for sm in SAMPLE_NAMES_INTERSECTION
    },
    **config["cdr_finder"],
}


# Avoid using github() due to Snakemake caching causing reruns.
module CDR_Finder:
    snakefile:
        "CDR-Finder/workflow/Snakefile"
    config:
        CDR_FINDER_CONFIG


use rule * from CDR_Finder as cdr_*


use rule calc_windows from CDR_Finder as cdr_calc_windows with:
    resources:
        mem=30,
        hrs=1,


use rule reorient_bed as reorient_cdr_bed with:
    input:
        bed=lambda wc: expand(rules.cdr_call_cdrs.output, sample=wc.sm),
        og_coords_key=rules.get_original_coords.output.og_coords_key,
    output:
        join(
            CDR_FINDER_OUTDIR,
            "bed",
            "{sm}_cdr_final.bed",
        ),
    log:
        join(CDR_FINDER_LOGDIR, "reorient_cdr_bed_{sm}.log"),
    params:
        legend_col_chrom_og="$1",
        legend_col_chrom_new='$4":"$5"-"$6',
        legend_col_chrom_len="$7",
        additional_cols="",


use rule reorient_bed as reorient_binned_methyl_bed with:
    input:
        bed=lambda wc: expand(
            rules.cdr_add_target_bed_coords_windows.output, sample=wc.sm
        ),
        og_coords_key=rules.get_original_coords.output.og_coords_key,
    output:
        temp(
            join(
                CDR_FINDER_OUTDIR,
                "bed",
                "{sm}_binned_freq_adj_final.bed",
            )
        ),
    log:
        join(CDR_FINDER_LOGDIR, "reorient_binned_methyl_bed_{sm}.log"),
    params:
        legend_col_chrom_og="$1",
        legend_col_chrom_new='$4":"$5"-"$6',
        legend_col_chrom_len="$7",
        additional_cols=", $4",


rule merge_cdr_beds:
    input:
        reorient_cdr_output=expand(
            rules.reorient_cdr_bed.output, sm=SAMPLE_NAMES_INTERSECTION
        ),
        reorient_methyl_cdr_output=expand(
            rules.reorient_binned_methyl_bed.output, sm=SAMPLE_NAMES_INTERSECTION
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
        cat {input.reorient_cdr_output} > {output.reorient_cdr_output}
        cat {input.reorient_methyl_cdr_output} > {output.reorient_methyl_cdr_output}
        """


rule cdr_finder_all:
    input:
        expand(rules.get_original_coords.output, sm=SAMPLE_NAMES_INTERSECTION),
        rules.merge_cdr_beds.output,
    default_target: True
