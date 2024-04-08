rule concat_asm:
    input:
        # Input directory per sample.
        sm_dir=os.path.join(config["concat_asm"]["input_dir"], "{sm}"),
    output:
        os.path.join(config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fasta"),
    # https://bioinf.shenwei.me/seqkit/usage/#rmdup
    params:
        assembly_fname_pattern="*.gz",
        by_seq="-s",
    resources:
        mem_mb=20_000,
    log:
        "logs/concat_asm_{sm}.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ zcat {input.sm_dir}/{params.assembly_fname_pattern} | \
        seqkit rmdup {params.by_seq};}} > {output} 2> {log}
        """


rule concat_asm_all:
    input:
        expand(rules.concat_asm.output, sm=SAMPLE_NAMES),
