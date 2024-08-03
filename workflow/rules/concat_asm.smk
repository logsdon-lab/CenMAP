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
        assembly_fname_pattern=".*\.(fa|fasta)",
    resources:
        mem=config["concat_asm"].get("mem", 4),
    log:
        "logs/concat_asm/concat_asm_{sm}.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ cat \
        <(find {input.sm_dir} -regextype posix-egrep -regex "{params.assembly_fname_pattern}\.gz" -size +0 -exec zcat {{}} + ) \
        <(find {input.sm_dir} -regextype posix-egrep -regex "{params.assembly_fname_pattern}" -size +0 -exec cat {{}} + ) | \
        seqkit rmdup;}} > {output.fa} 2> {log}
        samtools faidx {output.fa} 2>> {log}
        """


rule concat_asm_all:
    input:
        expand(rules.concat_asm.output, sm=SAMPLE_NAMES),
