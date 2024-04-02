# ex. NA12329 chrX    haplotype1      0000017 92036218        98212716        6176499 4617952 6176499 6176500
# ex. HG00171 chr22   h1tg000027l     1       26260313        1       4719177 4719177 48      4719177 4719178
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
        io_cols=" ".join(["ctg", "start", "end"]),
        grp_cols=" ".join(["ctg"]),
        sort_cols=" ".join(["ctg", "start"]),
    conda:
        "../env/py.yaml"
    log:
        "logs/make_{sm}_bed_files_for_plot.log",
    shell:
        # Only filter for sample to avoid malformed output ref cols in alr bed.
        # Also, filter starting position of 1 as likely only a fragment of ALR.
        """
        {{ cat {input.faidx} | \
        sed -e 's/_/\\t/g' -e 's/:/\\t/g' -e 's/-/\\t/g' -e 's/#/\\t/g' | \
        awk -v OFS="\\t" '{{
            if ($1 == "{wildcards.sm}") {{
                if ($3 ~ "h1" || $3 ~ "h2") {{
                    contig_name=$1"_"$2"_"$3"#"$4"-"$5
                    print contig_name, $6, $7
                }} else {{
                    contig_name=$1"_"$2"_"$3"-"$4
                    print contig_name, $5, $6
                }}
            }}
        }}' | \
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
        asm=rules.asm_rename_ctgs.output,
        reads=os.path.join(config["nuc_freq"]["hifi_reads_dir"], "{sm}", "{id}.bam"),
    output:
        alignment=temp(
            os.path.join(config["nuc_freq"]["output_dir"], "{sm}_{id}_hifi.bam")
        ),
    threads: config["nuc_freq"]["threads"]
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
        -j {threads} {input.asm} {input.reads} > {output.alignment} 2> {log}
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
    resources:
        sort_mem=4,
        mem_mb=20_000,
    threads: config["nuc_freq"]["threads"]
    conda:
        "../env/tools.yaml"
    log:
        "logs/filter_align_{sm}_{id}_hifi_reads_to_asm.log",
    shell:
        """
        {{ samtools view -b -F {params.samtools_view_flag} {input} | \
        samtools sort -m {resources.sort_mem}G -@ {threads} -o {output.alignment};}} 2> {log}
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
        mem_mb=90_000,
    params:
        ylim=100,
        height=4,
        # haplotype1 or h1
        hap_pattern=lambda wc: f"{wc.hap}|{str(wc.hap)[0]}{str(wc.hap)[-1]}",
    log:
        "logs/run_nucfreq_{sm}_{hap}.log",
    benchmark:
        "benchmarks/run_nucfreq_{sm}_{hap}.tsv"
    shell:
        """
        grep -P "{params.hap_pattern}" {input.alr_regions} > {output.alr_hap_regions}
        python {input.script} \
        -y {params.ylim} \
        {input.bam_file} \
        {output.plot} \
        --bed {output.alr_hap_regions} \
        --height {params.height} &> {log}
        """


# Then review plots manually.
