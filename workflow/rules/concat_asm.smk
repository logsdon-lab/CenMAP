rule concat_asm:
    input:
        # Input directory per sample.
        sm_dir=os.path.join(config["concat_asm"]["input_dir"], "{sm}"),
    output:
        fa=os.path.join(config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fa"),
        idx=os.path.join(
            config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fa.fai"
        ),
    # https://bioinf.shenwei.me/seqkit/usage/#rmdup
    params:
        assembly_fname_pattern=r".*\\.(fa|fasta)",
        assembly_fname_pattern_gz=r".*\\.(fa|fasta)\\.gz",
    resources:
        mem=config["concat_asm"].get("mem", "4GB"),
    log:
        "logs/concat_asm/concat_asm_{sm}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ cat \
        <(find {input.sm_dir} -regextype posix-egrep -regex "{params.assembly_fname_pattern_gz}" -size +0 -exec zcat {{}} + ) \
        <(find {input.sm_dir} -regextype posix-egrep -regex "{params.assembly_fname_pattern}" -size +0 -exec cat {{}} + ) | \
        seqkit rmdup;}} > {output.fa} 2> {log}
        samtools faidx {output.fa} 2>> {log}
        """


rule concat_asm_all:
    input:
        expand(rules.concat_asm.output, sm=SAMPLE_NAMES),
