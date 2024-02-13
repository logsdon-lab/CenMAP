rule concat_asm:
    input:
        # Input directory per sample.
        sm_dir=os.path.join(config["concat_asm"]["input_dir"], "{sm}"),
    output:
        os.path.join(
            config["concat_asm"]["output_dir"], "{sm}.vrk-ps-sseq.asm-comb-dedup.fasta"
        ),
    # https://bioinf.shenwei.me/seqkit/usage/#rmdup
    params:
        assembly_fname_pattern="*.gz",
        by_seq="-s",
    # Only allow not hap names ex. HG00171 (Not H00171_1)
    wildcard_constraints:
        sm="\\w+",
    log:
        "logs/concat_asm_{sm}.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ find {input.sm_dir} -name {params.assembly_fname_pattern} -exec zcat {{}} + | \
        seqkit rmdup {params.by_seq};}} > {output} 2> {log}
        """


rule concat_asm_all:
    input:
        expand(rules.concat_asm.output, sm=SAMPLE_NAMES),
