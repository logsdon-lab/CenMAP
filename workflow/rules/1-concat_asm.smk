include: "common.smk"
include: "utils.smk"


CONCAT_ASM_OUTDIR = join(OUTPUT_DIR, "1-concat_asm")
CONCAT_ASM_LOGDIR = join(LOG_DIR, "1-concat_asm")
CONCAT_ASM_BMKDIR = join(BMK_DIR, "1-concat_asm")


rule concat_asm:
    input:
        # Input directory per sample.
        sm_dir=join(config["concat_asm"]["input_dir"], "{sm}"),
    output:
        fa=join(CONCAT_ASM_OUTDIR, "{sm}-asm-comb-dedup.fa"),
        idx=join(CONCAT_ASM_OUTDIR, "{sm}-asm-comb-dedup.fa.fai"),
    # Remove duplicate contigs. https://bioinf.shenwei.me/seqkit/usage/#rmdup
    # Remove descriptions. https://bioinf.shenwei.me/seqkit/usage/#replace
    params:
        assembly_fname_pattern=r".*\\.(fa|fasta)",
        assembly_fname_pattern_gz=r".*\\.(fa|fasta)\\.gz",
    resources:
        mem=config["concat_asm"].get("mem", "4GB"),
    log:
        join(CONCAT_ASM_LOGDIR, "concat_asm_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ cat \
        <(find {input.sm_dir}/ -regextype posix-egrep -regex "{params.assembly_fname_pattern_gz}" -size +0 -exec zcat {{}} + ) \
        <(find {input.sm_dir}/ -regextype posix-egrep -regex "{params.assembly_fname_pattern}" -size +0 -exec cat {{}} + ) | \
        seqkit rmdup | \
        seqkit replace -p "\\s.+" ;}} > {output.fa} 2> {log}
        samtools faidx {output.fa} 2>> {log}
        """


rule concat_asm_all:
    input:
        expand(rules.concat_asm.output, sm=SAMPLE_NAMES),
    default_target: True
