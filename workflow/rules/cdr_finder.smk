include: "common.smk"


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


rule merge_methyl_bam_to_fq:
    input:
        os.path.join(config["cdr_finder"]["input_bam_dir"], "{sm}"),
    output:
        temp(
            os.path.join(
                config["cdr_finder"]["output_dir"], "aln", "{sm}_methyl.fq.gz"
            )
        ),
    params:
        file_pattern=config["cdr_finder"]["file_pattern"],
    resources:
        mem="16G",
        sort_mem="8G",
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
        samtools bam2fq -T "*" -@ {threads} - | \
        bgzip ;}} > {output} 2> {log}
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
        kmer_cnts=rules.get_kmer_cnts.output.kmer_cnts_list,
    output:
        bam=os.path.join(config["cdr_finder"]["output_dir"], "aln", "{sm}.bam"),
    params:
        preset="map-ont",
        samtools_view_flag=2038,
        min_peak_dp_aln_score=4000,
        split_idx_num_base="10G",
    threads: config["cdr_finder"]["aln_threads"]
    resources:
        mem=config["cdr_finder"]["aln_mem"],
    conda:
        "../envs/winnowmap.yaml"
    log:
        "logs/cdr_finder/align_methyl_bam_to_asm_{sm}.log",
    benchmark:
        "benchmarks/cdr_finder/align_methyl_bam_to_asm_{sm}.tsv"
    shell:
        """
        {{ winnowmap -W {input.kmer_cnts} \
        -y --eqx \
        -ax {params.preset} \
        -s {params.min_peak_dp_aln_score} \
        -t {threads} \
        -I {params.split_idx_num_base} \
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
        join -1 3 -2 1 \
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
            print $1, adj_start, adj_end, ctg, start, end
        }}' > {output.og_coords_key}
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


module CDR_Finder:
    snakefile:
        github(
            "koisland/CDR-Finder",
            path="workflow/Snakefile",
            branch="feature/add-test-data",
        )
    config:
        CDR_FINDER_CONFIG


use rule * from CDR_Finder as cdr_*


use rule calc_windows from CDR_Finder as cdr_calc_windows with:
    resources:
        mem=30,
        hrs=1,


rule cdr_finder_only:
    input:
        expand(rules.cdr_all.input, sm=SAMPLE_NAMES_INTERSECTION),
