

rule convert_reads_to_fq:
    input:
        reads=ancient(
            os.path.join(config["nuc_freq"]["hifi_reads_dir"], "{sm}", "{id}")
        ),
    output:
        temp(os.path.join(config["nuc_freq"]["output_dir"], "{sm}_{id}.fq")),
    conda:
        "../env/tools.yaml"
    log:
        "logs/convert_{sm}_{id}_to_fq.log",
    shell:
        """
        if [[ "{wildcards.id}" =~ .*\.bam$ ]]; then
            samtools bam2fq {input.reads} > {output} 2> {log}
        elif [[ "{wildcards.id}" =~ .*\.gz$ ]]; then
            zcat {input.reads} > {output} 2> {log}
        else
            cp {input.reads} {output} 2> {log}
        fi
        """


rule align_reads_to_asm:
    input:
        asm=ancient(
            os.path.join(
                config["concat_asm"]["output_dir"],
                "{sm}",
                "{sm}_regions.renamed.fa",
            )
        ),
        reads=ancient(rules.convert_reads_to_fq.output),
    output:
        temp(os.path.join(config["nuc_freq"]["output_dir"], "{sm}_{id}_hifi.bam")),
    threads: config["nuc_freq"]["threads_aln"]
    resources:
        mem_mb=120_000,
    params:
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
        pbmm2 align \
        --log-level {params.aln_log_level} \
        --preset {params.aln_preset} \
        --min-length {params.aln_min_length} \
        -j {threads} {input.asm} {input.reads} > {output} 2>> {log}
        """


# Get error when trying to pipe ^ to samtools view. No header. Separate step works.
rule filter_align_reads_to_asm:
    input:
        ancient(rules.align_reads_to_asm.output),
    output:
        temp(os.path.join(config["nuc_freq"]["output_dir"], "{sm}_{id}_hifi_view.bam")),
    params:
        # https://broadinstitute.github.io/picard/explain-flags.html
        samtools_view_flag=config["nuc_freq"]["samtools_view_flag"],
    resources:
        sort_mem=4,
        mem_mb=20_000,
    threads: config["nuc_freq"]["threads_aln"]
    conda:
        "../env/tools.yaml"
    log:
        "logs/filter_align_{sm}_{id}_hifi_reads_to_asm.log",
    shell:
        """
        {{ samtools view -b -F {params.samtools_view_flag} {input} | \
        samtools sort -m {resources.sort_mem}G -@ {threads} -o {output};}} 2> {log}
        """


def get_aln_to_asm(wc) -> list[str]:
    alns = expand(
        rules.filter_align_reads_to_asm.output,
        sm=[wc.sm],
        id=SAMPLE_FLOWCELL_IDS[str(wc.sm)],
    )
    if not alns:
        raise FileNotFoundError(
            f"Subdirectory {wc.sm} in {config['nuc_freq']['hifi_reads_dir']} is missing or contains no alignment files."
        )
    return ancient(alns)


rule merge_hifi_read_asm_alignments:
    input:
        get_aln_to_asm,
    output:
        alignment=os.path.join(config["nuc_freq"]["output_dir"], "{sm}_hifi.bam"),
        alignment_idx=os.path.join(
            config["nuc_freq"]["output_dir"], "{sm}_hifi.bam.bai"
        ),
    threads: config["nuc_freq"]["threads_aln"]
    resources:
        mem_mb=10_000,
    conda:
        "../env/tools.yaml"
    log:
        "logs/merge_{sm}_hifi_read_asm_alignments.log",
    benchmark:
        "benchmarks/merge_{sm}_hifi_read_asm_alignments.tsv"
    shell:
        """
        samtools merge -@ {threads} {output.alignment} {input} 2> {log}
        samtools index {output.alignment} 2> {log}
        """


rule check_asm_nucfreq:
    input:
        bam_file=ancient(rules.merge_hifi_read_asm_alignments.output.alignment),
        alr_regions=rules.make_new_cens_bed_file.output.alr_bed,
        config=config["nuc_freq"]["config_nucfreq"],
        ignore_regions=config["nuc_freq"]["ignore_regions"],
    output:
        plot_dir=directory(os.path.join(config["nuc_freq"]["output_dir"], "{sm}")),
        misassemblies=os.path.join(
            config["nuc_freq"]["output_dir"],
            "{sm}_cen_misassemblies.bed",
        ),
        asm_status=os.path.join(
            config["nuc_freq"]["output_dir"],
            "{sm}_cen_status.bed",
        ),
    threads: config["nuc_freq"]["processes_nucfreq"]
    conda:
        "../env/nucfreq.yaml"
    resources:
        mem_mb=50_000,
    log:
        "logs/run_nucfreq_{sm}.log",
    benchmark:
        "benchmarks/run_nucfreq_{sm}.tsv"
    shell:
        """
        nucfreq \
        -i {input.bam_file} \
        -b {input.alr_regions} \
        -d {output.plot_dir} \
        -o {output.misassemblies} \
        -t {threads} \
        -p {threads} \
        -s {output.asm_status} \
        -c {input.config} \
        --ignore_regions {input.ignore_regions} &> {log}
        """


rule nuc_freq_only:
    input:
        expand(rules.merge_hifi_read_asm_alignments.output, sm=SAMPLE_NAMES),
        expand(rules.check_asm_nucfreq.output, sm=SAMPLE_NAMES),
