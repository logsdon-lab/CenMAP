# HG00171_chr4_haplotype1-0000002:1892469-12648706
# |haplotype1-0000002|1892469|12648706|chr4|
# |haplotype1-0000002|1892469|12648706|12648706-1892469|chr4|
rule make_bed_files_for_plot:
    input:
        script="workflow/scripts/filter_cen_ctgs.py",
        faidx=expand(
            rules.new_cens_index_renamed_ctgs.output,
            ort=ORIENTATION,
            sm=SAMPLE_NAMES,
        ),
    output:
        os.path.join(config["nuc_freq"]["output_dir"], "{sm}_ALR_regions.500kp.bed"),
    params:
        io_cols=" ".join(["ctg", "start", "end", "chr"]),
        grp_cols=" ".join(["ctg", "chr"]),
        sort_cols=" ".join(["ctg", "start"]),
    conda:
        "../env/py.yaml"
    log:
        "logs/make_{sm}_bed_files_for_plot.log",
    shell:
        """
        {{ cat {input.faidx} | \
        sed -e 's/_/\\t/g' -e 's/:/\\t/g' -e 's/-/\\t/g' | \
        awk -v OFS="\\t" '{{print $3"-"$4, $5, $6, $2}}' | \
        sort | \
        uniq | \
        python {input.script} bedminmax \
            -ci {params.io_cols} \
            -co {params.io_cols} \
            -g {params.grp_cols} \
            -s {params.sort_cols} | \
        awk -v OFS="\\t" '{{print $1, $2, $3, $3-$2, $4}}';}} > {output} 2> {log}
        """


# rule convert_hifi_reads_to_fq:
#     input:
#         lambda wc: expand(
#             os.path.join(reads_dir, wc.sm, "{id}.bam"),
#             id=SAMPLE_FLOWCELL_IDS[wc.sm]
#         ),
#     output:
#         os.path.join(config["nuc_freq"]["output_dir"], "{sm}_{id}_hifi.fq"),
#     conda:
#         "../env/tools.yaml"
#     log:
#         "logs/convert_{sm}_{id}_hifi_reads_to_fq.log",
#     shell:
#         """
#         samtools bam2fq {input} > {output} 2> {log}
#         """


rule align_reads_to_asm:
    input:
        asm=rules.concat_asm.output,
        reads=os.path.join(config["nuc_freq"]["hifi_reads_dir"], "{sm}", "{id}.bam"),
    output:
        alignments=temp(
            os.path.join(config["nuc_freq"]["output_dir"], "{sm}_{id}_hifi.bam")
        ),
    threads: config["nuc_freq"]["threads"]
    resources:
        sort_mem=4,
    params:
        # https://broadinstitute.github.io/picard/explain-flags.html
        samtools_view_flag=config["nuc_freq"]["samtools_view_flag"],
        aln_log_level="DEBUG",
        aln_preset="SUBREAD",
        aln_min_length=5000,
    conda:
        "../env/tools.yaml"
    log:
        "logs/align_{sm}_{id}_hifi_reads_to_asm.log",
    benchmark:
        "benchmarks/align_{sm}_{id}_hifi_reads_to_asm.tsv"
    shell:
        """
        {{ pbmm2 align \
        --log-level {params.aln_log_level} \
        --preset {params.aln_preset} \
        --min-length {params.aln_min_length} \
        -j {threads} {input.asm} {input.reads} | \
        samtools view -F {params.samtools_view_flag} -u - | \
        samtools sort -m {resources.sort_mem}G -@ {threads} -;}} > {output} 2> {log}
        """


rule merge_hifi_read_asm_alignments:
    input:
        lambda wc: expand(
            rules.align_reads_to_asm.output,
            sm=[wc.sm],
            id=SAMPLE_FLOWCELL_IDS[str(wc.sm)],
        ),
    output:
        alignment=os.path.join(config["nuc_freq"]["output_dir"], "{sm}_hifi.bam"),
    conda:
        "../env/tools.yaml"
    log:
        "logs/merge_{sm}_hifi_read_asm_alignments.log",
    benchmark:
        "benchmarks/merge_{sm}_hifi_read_asm_alignments.tsv"
    shell:
        """
        samtools merge -@ {threads} {output} {input} 2> {log}
        """


# Merge now or pass each one to vvv?
rule gen_nucfreq_plot:
    input:
        script="workflow/scripts/NucPlot.py",
        bam_file=rules.merge_hifi_read_asm_alignments.output,
        alr_regions=rules.make_bed_files_for_plot.output,
    output:
        plot=os.path.join(
            config["nuc_freq"]["output_dir"],
            "{sm}_hifi_cens.png",
        ),
    conda:
        "../env/py.yaml"
    params:
        ylim=100,
        height=4,
    log:
        "logs/run_nucfreq_{sm}.log",
    benchmark:
        "benchmarks/run_nucfreq_{sm}.tsv"
    shell:
        """
        python {input.script} \
        -y {params.ylim} \
        {input.bam_file} \
        {output} \
        --bed {input.alr_regions} \
        --height {params.height} &> {log}
        """


# Then review plots manually.
