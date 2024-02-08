rule concat_asm:
    input:
        expand(
            os.path.join(
                config["concat_asm"]["input_dir"],
                "{{sm}}",
                    "{{sm}}.vrk-ps-sseq.{typ}.fasta.gz",
                ),
                typ=[
                typ if typ == "contaminants" else f"asm-{typ}"
                for typ in config["concat_asm"]["types"]
            ],
        ),
    output:
        os.path.join(
            config["concat_asm"]["output_dir"], "{sm}.vrk-ps-sseq.asm-comb-dedup.fasta"
        ),
    # https://bioinf.shenwei.me/seqkit/usage/#rmdup
    params:
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
        {{ zcat {input} | seqkit rmdup {params.by_seq};}} > {output} 2> {log}
        """


rule concat_asm_all:
    input:
        expand(rules.concat_asm.output, sm=SAMPLE_NAMES),
