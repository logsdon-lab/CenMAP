include: "common.smk"
include: "utils.smk"


# Get sample names with subdirs in cdr_finder.input_bam_dir
SAMPLE_NAMES_BAM = set(
    glob_wildcards(os.path.join(config["cdr_finder"]["input_bam_dir"], "{sm}")).sm
)
SAMPLE_NAMES_INTERSECTION = SAMPLE_NAMES_BAM.intersection(set(SAMPLE_NAMES))
print(
    f"Using {len(SAMPLE_NAMES_INTERSECTION)} out of {len(SAMPLE_NAMES)} samples "
    f"that have subdirs in {config['cdr_finder']['input_bam_dir']} for CDR-Finder.",
    file=sys.stderr,
)
ALIGNER = config["cdr_finder"].get("aligner", "minimap2")
ADDED_ALIGNER_OPTS = config["cdr_finder"].get("aligner_added_opts", "")
ALL_ALIGNER_SETTINGS = {
    "minimap2": {
        "preset": "lr:hqae",
        "split_idx_num_base": "-I 10G",
        "min_peak_dp_aln_score": "-s 4000",
        "additional": ADDED_ALIGNER_OPTS,
    },
    # https://www.biorxiv.org/content/10.1101/2024.11.01.621587v1
    "winnowmap": {
        "preset": "map-ont",
        "split_idx_num_base": "-I 10G",
        "min_peak_dp_aln_score": "-s 4000",
        "additional": "-W {input.kmer_cnts} " + ADDED_ALIGNER_OPTS,
    },
}
try:
    ALIGNER_SETTINGS = ALL_ALIGNER_SETTINGS[ALIGNER]
except KeyError:
    raise ValueError(
        f"Invalid aligner option ({ALIGNER}) for CDR-Finder. Choose: {tuple(ALL_ALIGNER_SETTINGS.keys())} "
    )


wildcard_constraints:
    sm="|".join(SAMPLE_NAMES_INTERSECTION),


rule merge_methyl_bam_to_fq:
    input:
        os.path.join(config["cdr_finder"]["input_bam_dir"], "{sm}"),
    output:
        temp(
            os.path.join(
                config["cdr_finder"]["output_dir"], "aln", "{sm}_methyl.fq"
            )
        ),
    params:
        file_pattern=config["cdr_finder"]["file_pattern"],
    resources:
        mem="16G",
    threads: config["cdr_finder"]["aln_threads"]
    log:
        "logs/cdr_finder/merge_methyl_bam_{sm}.log",
    benchmark:
        "benchmarks/cdr_finder/merge_methyl_bam_{sm}.tsv"
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ samtools merge -@ {threads} - $(find {input} -regex {params.file_pattern}) | \
        samtools bam2fq -T "*" -@ {threads} - ;}} > {output} 2> {log}
        """


rule get_kmer_cnts:
    input:
        fa=os.path.join(config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fa"),
    output:
        kmer_cnts=directory(
            os.path.join(config["cdr_finder"]["output_dir"], "aln", "{sm}")
        ),
        kmer_cnts_list=os.path.join(
            config["cdr_finder"]["output_dir"], "aln", "{sm}_repetitive.txt"
        ),
    params:
        kmer_size=15,
        gt_perc_thr=0.9998,
    threads: config["cdr_finder"]["aln_threads"]
    resources:
        mem="16G",
    log:
        "logs/cdr_finder/get_kmer_cnts_{sm}.log",
    benchmark:
        "benchmarks/cdr_finder/get_kmer_cnts_{sm}.tsv"
    conda:
        "../envs/winnowmap.yaml"
    shell:
        """
        meryl count threads={threads} k={params.kmer_size} output {output.kmer_cnts} {input.fa} 2> {log}
        meryl print greater-than distinct={params.gt_perc_thr} {output.kmer_cnts} > {output.kmer_cnts_list} 2>> {log}
        """


rule align_methyl_bam_to_asm:
    input:
        ref=os.path.join(config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fa"),
        query=rules.merge_methyl_bam_to_fq.output,
        kmer_cnts=(
            rules.get_kmer_cnts.output.kmer_cnts_list if ALIGNER == "winnowmap" else []
        ),
    output:
        bam=os.path.join(config["cdr_finder"]["output_dir"], "aln", "{sm}.bam"),
    params:
        aligner=ALIGNER,
        aligner_added_opts=ADDED_ALIGNER_OPTS,
        preset=ALIGNER_SETTINGS["preset"],
        samtools_view_flag=2038,
        min_peak_dp_aln_score=ALIGNER_SETTINGS["min_peak_dp_aln_score"],
        split_idx_num_base=ALIGNER_SETTINGS["split_idx_num_base"],
    threads: config["cdr_finder"]["aln_threads"]
    resources:
        mem=config["cdr_finder"]["aln_mem"],
    conda:
        f"../envs/{ALIGNER}.yaml"
    log:
        "logs/cdr_finder/align_methyl_bam_to_asm_{sm}.log",
    benchmark:
        "benchmarks/cdr_finder/align_methyl_bam_to_asm_{sm}.tsv"
    shell:
        """
        {{ {params.aligner} \
        -y --eqx \
        -ax {params.preset} \
        -t {threads} \
        {params.min_peak_dp_aln_score} \
        {params.split_idx_num_base} \
        {params.aligner_added_opts} \
        {input.ref} {input.query} | \
        samtools view -u -F {params.samtools_view_flag} - | \
        samtools sort -o {output.bam} ;}} 2> {log}
        """


rule get_original_coords:
    input:
        bed=os.path.join(
            config["repeatmasker"]["output_dir"],
            "bed",
            "{sm}_complete_correct_ALR_regions.bed",
        ),
        asm_faidx=os.path.join(
            config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fa.fai"
        ),
    output:
        og_coords=os.path.join(
            config["cdr_finder"]["output_dir"],
            "bed",
            "{sm}_complete_correct_ALR_regions.og_coords.bed",
        ),
        og_coords_key=os.path.join(
            config["cdr_finder"]["output_dir"],
            "bed",
            "{sm}_complete_correct_ALR_regions.og_coords.key.bed",
        ),
    log:
        "logs/cdr_finder/get_original_coords_{sm}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ join -1 3 -2 1 \
            <(sed 's/_/\\t/g' {input.bed} | sort -k 3 | cut -f 1,2,3,4,5) \
            <(cut -f 1,2 {input.asm_faidx} | sort -k 1) | \
        awk -v OFS="\\t" '{{
            ctg=$2"_"$3"_"$1
            is_reversed=($3 ~ "rc")
            start=$4; end=$5;
            adj_start=$4; adj_end=$5;
            if (is_reversed) {{
                adj_start=$6-end
                adj_end=$6-start
            }}
            print $1, adj_start, adj_end > "{output.og_coords}"
            print $1, adj_start, adj_end, ctg, start, end, $6
        }}';}} > {output.og_coords_key} 2> {log}
        """


# Pass CDR config here.
CDR_FINDER_CONFIG = {
    **config["cdr_finder"],
    "samples": {
        sm: {
            "fasta": os.path.join(
                config["concat_asm"]["output_dir"], f"{sm}-asm-comb-dedup.fa"
            ),
            "regions": expand(rules.get_original_coords.output.og_coords, sm=sm),
            "bam": expand(rules.align_methyl_bam_to_asm.output, sm=sm),
        }
        for sm in SAMPLE_NAMES_INTERSECTION
    },
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
        os.path.join(
            config["cdr_finder"]["output_dir"],
            "bed",
            "{sm}_cdr_final.bed",
        ),
    log:
        "logs/cdr_finder/reorient_cdr_bed_{sm}.log",
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
            os.path.join(
                config["cdr_finder"]["output_dir"],
                "bed",
                "{sm}_binned_freq_adj_final.bed",
            )
        ),
    log:
        "logs/cdr_finder/reorient_binned_methyl_bed_{sm}.log",
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
        reorient_cdr_output=os.path.join(
            config["cdr_finder"]["output_dir"],
            "bed",
            "all_cdrs.bed",
        ),
        reorient_methyl_cdr_output=os.path.join(
            config["cdr_finder"]["output_dir"],
            "bed",
            "all_binned_freq.bed",
        ),
    shell:
        """
        cat {input.reorient_cdr_output} > {output.reorient_cdr_output}
        cat {input.reorient_methyl_cdr_output} > {output.reorient_methyl_cdr_output}
        """


rule cdr_finder_only:
    input:
        expand(rules.get_original_coords.output, sm=SAMPLE_NAMES_INTERSECTION),
        rules.merge_cdr_beds.output,
