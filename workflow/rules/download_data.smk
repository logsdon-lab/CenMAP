include: "common.smk"


ASMS_INPUT_DIR = config["concat_asm"]["input_dir"]
READS_INPUT_DIR = config["nucflag"]["hifi_reads_dir"]
ASM_URIS = {"HG002": "s3://human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz"}
HIFI_URIS = {"HG002": "s3://human-pangenomics/T2T/scratch/HG002/sequencing/hifirevio/"}


rule download_assembly:
    output:
        directory(os.path.join(ASMS_INPUT_DIR, "{sm}")),
    params:
        uri=lambda wc: ASM_URIS[str(wc.sm)],
    conda:
        "../envs/nucflag_bmk.yaml"
    log:
        "logs/download_data/download_assembly_{sm}.log",
    shell:
        """
        aws s3 --no-sign-request cp {params.uri} {output} 2> {log}
        """


rule download_hifi:
    output:
        directory(os.path.join(READS_INPUT_DIR, "{sm}")),
    threads: 20
    params:
        uri=lambda wc: HIFI_URIS[str(wc.sm)],
    conda:
        "../envs/nucflag_bmk.yaml"
    log:
        "logs/download_data/download_hifi_{sm}.log",
    shell:
        """
        aws s3 --no-sign-request sync {params.uri} {output} 2> {log}
        """


checkpoint download_data_all:
    input:
        expand(rules.download_assembly.output, sm=SAMPLE_NAMES),
        expand(rules.download_hifi.output, sm=SAMPLE_NAMES),
    output:
        touch("download_data.done"),


rule download_bmk_data_all:
    input:
        rules.download_data_all.output,
    default_target: True
