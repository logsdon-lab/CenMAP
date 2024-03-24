# HG00171_chr4_haplotype1-0000002:1892469-12648706
# |haplotype1-0000002|1892469|12648706|chr4|
# |haplotype1-0000002|1892469|12648706|12648706-1892469|chr4|
rule make_bed_files_for_plot:
    input:
        script="workflow/scripts/filter_cen_ctgs.py",
        faidx=expand(
            rules.merge_alr_regions_by_chr.output.idx,
            ort=ORIENTATION,
            sm=SAMPLE_NAMES,
            chr=CHROMOSOMES,
        ),
    output:
        tmp_fmt_alr_bed=temp(
            os.path.join(
                config["nuc_freq"]["output_dir"], "fmt_{sm}_ALR_regions.500kbp.bed"
            )
        ),
        alr_bed=os.path.join(
            config["nuc_freq"]["output_dir"], "{sm}_ALR_regions.500kbp.bed"
        ),
        correct_alr_bed=os.path.join(
            config["nuc_freq"]["output_dir"], "{sm}_correct_ALR_regions.500kbp.bed"
        ),
    params:
        io_cols=" ".join(["ctg", "start", "end", "chr"]),
        grp_cols=" ".join(["ctg", "chr"]),
        sort_cols=" ".join(["ctg", "start"]),
    conda:
        "../env/py.yaml"
    log:
        "logs/make_{sm}_bed_files_for_plot.log",
    message:
        f"""
        Check (sm)_correct_ALR_regions.500kbp.bed for copies of ALR regions per sample.
        Add a column labeling regions as "good" or "misassembled".
        Then update repeatmasker.correct_asm in config/config.yaml.
        An automated approach to missassembly identification is a WIP.
        """
    shell:
        # Only filter for sample to avoid malformed output ref cols in alr bed.
        # Also, filter starting position of 1 as likely only a fragment of ALR.
        """
        {{ cat {input.faidx} | \
        sed -e 's/_/\\t/g' -e 's/:/\\t/g' -e 's/-/\\t/g' | \
        awk -v OFS="\\t" '{{if ($1 == "{wildcards.sm}") {{print $3"-"$4, $5, $6, $2}}}}' | \
        sort | \
        uniq;}} > {output.tmp_fmt_alr_bed} 2> {log}

        {{ python {input.script} bedminmax \
            -i {output.tmp_fmt_alr_bed} \
            -ci {params.io_cols} \
            -co {params.io_cols} \
            -g {params.grp_cols} \
            -s {params.sort_cols} | \
        awk -v OFS="\\t" '{{ if ($2 != 1) {{print $1, $2, $3, $3-$2, $4}}}}';}} > {output.alr_bed} 2>> {log}

        # Make copy.
        cp {output.alr_bed} {output.correct_alr_bed}
        """


rule align_reads_to_asm:
    input:
        asm=rules.concat_asm.output,
        reads=os.path.join(config["nuc_freq"]["hifi_reads_dir"], "{sm}", "{id}.bam"),
    output:
        alignment=temp(
            os.path.join(config["nuc_freq"]["output_dir"], "{sm}_{id}_hifi.bam")
        ),
    threads: config["nuc_freq"]["threads"]
    resources:
        mem_mb=180_000,
        sort_mem=4,
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
        -j {threads} {input.asm} {input.reads} \
        --sort \
        --sort-memory {resources.sort_mem}G \
        --sort-threads {threads} > {output.alignment} 2> {log}
        """


# Get error when trying to pipe ^ to samtools view. No header. Separate step works.
rule filter_align_reads_to_asm:
    input:
        rules.align_reads_to_asm.output,
    output:
        alignment=temp(
            os.path.join(config["nuc_freq"]["output_dir"], "{sm}_{id}_hifi_view.bam")
        ),
    params:
        # https://broadinstitute.github.io/picard/explain-flags.html
        samtools_view_flag=config["nuc_freq"]["samtools_view_flag"],
    conda:
        "../env/tools.yaml"
    log:
        "logs/filter_align_{sm}_{id}_hifi_reads_to_asm.log",
    shell:
        """
        samtools view -bo {output.alignment} -F {params.samtools_view_flag} {input} 2> {log}
        """


rule merge_hifi_read_asm_alignments:
    input:
        lambda wc: expand(
            rules.filter_align_reads_to_asm.output.alignment,
            sm=[wc.sm],
            id=SAMPLE_FLOWCELL_IDS[str(wc.sm)],
        ),
    output:
        alignment=os.path.join(config["nuc_freq"]["output_dir"], "{sm}_hifi.bam"),
        alignment_idx=os.path.join(
            config["nuc_freq"]["output_dir"], "{sm}_hifi.bam.bai"
        ),
    threads: config["nuc_freq"]["threads"]
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


rule gen_nucfreq_plot:
    input:
        script="workflow/scripts/NucFreq/NucPlot.py",
        bam_file=rules.merge_hifi_read_asm_alignments.output.alignment,
        alr_regions=rules.make_bed_files_for_plot.output.alr_bed,
    output:
        alr_hap_regions=temp(
            os.path.join(
                config["nuc_freq"]["output_dir"],
                "{sm}_{hap}_ALR_regions.500kbp.bed",
            )
        ),
        plot=os.path.join(
            config["nuc_freq"]["output_dir"],
            "{sm}_{hap}_hifi_cens.png",
        ),
    conda:
        "../env/pysam.yaml"
    resources:
        mem_mb=60_000,
    params:
        ylim=100,
        height=4,
    log:
        "logs/run_nucfreq_{sm}_{hap}.log",
    benchmark:
        "benchmarks/run_nucfreq_{sm}_{hap}.tsv"
    shell:
        """
        grep "{wildcards.hap}" {input.alr_regions} > {output.alr_hap_regions}
        python {input.script} \
        -y {params.ylim} \
        {input.bam_file} \
        {output.plot} \
        --bed {output.alr_hap_regions} \
        --height {params.height} &> {log}
        """


# Then review plots manually.
