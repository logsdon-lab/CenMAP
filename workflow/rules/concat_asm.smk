CONCAT_ASM_INPUT_DIR = config["concat_asm"]["input_dir"]
CONCAT_ASM_OUTPUT_DIR = config["concat_asm"]["output_dir"]


rule concat_asm:
    input:
        expand(
            os.path.join(
                CONCAT_ASM_INPUT_DIR, "{{sm}}", "{{sm}}.vrk-ps-sseq.asm-{typ}.fasta.gz"
            ),
            typ=config["concat_asm"]["types"],
        ),
    output:
        os.path.join(CONCAT_ASM_OUTPUT_DIR, "{sm}.vrk-ps-sseq.asm-comb-dedup.fa.gz"),
    # https://bioinf.shenwei.me/seqkit/usage/#rmdup
    params:
        by_seq="-s",
    # Only allow not hap names ex. HG00171 (Not H00171_1)
    wildcard_constraints:
        sm="\w+",
    log:
        "logs/concat_asm_{sm}.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ zcat {input} | seqkit rmdup {params.by_seq} | gzip;}} > {output} 2> {log}
        """


rule concat_asm_all:
    input:
        expand(rules.concat_asm.output, sm=SAMPLE_NAMES),
