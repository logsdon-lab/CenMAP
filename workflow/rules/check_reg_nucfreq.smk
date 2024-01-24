# HG00171_chr4_haplotype1-0000002:1892469-12648706
# |haplotype1-0000002|1892469|12648706|chr4|
# |haplotype1-0000002|1892469|12648706|12648706-1892469|chr4|
rule make_bed_files_for_plot:
    input:
        script="workflow/scripts/filter_cen_ctgs.py",
        faidx=lambda wc: expand(
            rules.new_cens_index_renamed_ctgs.output, ort=ORIENTATION, sm=[wc.sm]
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


def get_hifi_read_files(wc) -> list[str]:
    """
    Get hifi reads by sample automatically from hifi_reads_dir.
    Expects {hifi_reads_dir}/{sample}/*.bam
    """
    path_pattern = os.path.join(
        config["nuc_freq"]["hifi_reads_dir"], wc.sm, "{mdata_id}.bam"
    )
    reads_run_mdata_id = glob_wildcards(path_pattern)
    return expand(path_pattern, mdata_id=reads_run_mdata_id.mdata_id)


rule convert_hifi_reads_to_fq:
    input:
        get_hifi_read_files,
    output:
        os.path.join(config["nuc_freq"]["output_dir"], "{sm}_hifi.fq"),
    conda:
        "../env/tools.yaml"
    log:
        "logs/convert_{sm}_hifi_reads_to_fq.log",
    shell:
        """
        for file in {input}; do
            samtools bam2fq $file >> {output} 2> {log}
        done
        """


rule align_reads_to_asm:
    input:
        asm=lambda wc: expand(rules.concat_asm.output, sm=[str(wc.sm).split("_")[0]]),
        reads=rules.convert_hifi_reads_to_fq.output,
    output:
        alignments=os.path.join(config["nuc_freq"]["output_dir"], "{sm}_hifi.bam"),
    threads: 20
    conda:
        "../env/tools.yaml"
    log:
        "logs/align_{sm}_hifi_reads_to_ref.log",
    shell:
        """
        {{ minimap2 -ax map-pb -t {threads} {input.asm} {input.reads} | \
        samtools view;}} > {output} 2> {log}
        """


# Merge now or pass each one to vvv?
rule gen_nucfreq_plot:
    input:
        script="workflow/scripts/NucPlot.py",
        bam_file=rules.align_reads_to_asm.output,
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


rule nuc_freq_all:
    input:
        expand(rules.gen_nucfreq_plot.output, sm=SAMPLES_DF.index),
