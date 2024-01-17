import os


rule concat_asm:
    input:
        expand(os.path.join(), types=config["concat_asm"]["types"]),
    output:
        "",
    shell:
        """
        zcat \
        data/HG00171/HG00171.vrk-ps-sseq.asm-hap1.fasta.gz \
        data/HG00171/HG00171.vrk-ps-sseq.asm-hap2.fasta.gz > HG000171.vrk-ps-sseq.asm-comb.fasta.gz
        """
