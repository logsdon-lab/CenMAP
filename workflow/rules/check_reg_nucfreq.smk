
rule align_reads_to_asm:
    input:
        asm=ancient(
            os.path.join(
                config["concat_asm"]["output_dir"],
                "{sm}",
                "{sm}_regions.renamed.fa",
            )
        ),
        reads=ancient(
            os.path.join(
                config["nuc_freq"]["hifi_reads_dir"],
                "{sm}",
                f"{{id}}.{config['nuc_freq']['reads_ext']}",
            )
        ),
    output:
        temp(os.path.join(config["nuc_freq"]["output_dir"], "{sm}_{id}_hifi.bam")),
    threads: config["nuc_freq"]["threads_aln"]
    resources:
        mem_mb=config["nuc_freq"]["mem_mb_aln"],
        sort_mem=4,
    params:
        aln_log_level="DEBUG",
        aln_preset="SUBREAD",
        aln_min_length=5000,
        tmp_dir=config["nuc_freq"].get("tmp_dir", os.environ.get("TMPDIR", "/tmp")),
        samtools_view_flag=config["nuc_freq"]["samtools_view_flag"],
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
        samtools sort -T {params.tmp_dir} -m {resources.sort_mem}G -@ {threads} - ;}} > {output} 2>> {log}
        """


def get_aln_to_asm(wc) -> list[str]:
    alns = expand(
        rules.align_reads_to_asm.output,
        sm=[wc.sm],
        id=SAMPLE_FLOWCELL_IDS[str(wc.sm)],
    )
    if not alns:
        raise FileNotFoundError(
            f"Directory {config['nuc_freq']['hifi_reads_dir']}/{wc.sm} is missing or contains no alignment files with extension {config['nuc_freq']['reads_ext']}."
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
        mem_mb=config["nuc_freq"]["mem_mb_nucfreq"],
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
        expand(rules.check_asm_nucfreq.output, sm=SAMPLE_NAMES),
