include: "common.smk"


READS_SAMPLE_SHEET = "config/b2_formatted.tsv"
READS_INPUT_DIR = config["nucflag"]["hifi_reads_fofn_dir"]
ASMS_INPUT_DIR = config["concat_asm"]["input_dir"]
SAMPLE_READ_URIS = defaultdict(list)


wildcard_constraints:
    sm="|".join(SAMPLE_NAMES),


with open(READS_SAMPLE_SHEET) as fh:
    for line in fh.readlines():
        sm, s3_uri = line.strip().split("\t")
        SAMPLE_READ_URIS[sm].append(s3_uri)


rule download_assemblies:
    output:
        directory(os.path.join(ASMS_INPUT_DIR, "{sm}")),
    params:
        verkko_thic_uri="s3://human-pangenomics/submissions/0624338D-1A6F-4E29-A276-D2C247FE0558--verkko-v2.1_intermediate_asms/{sm}/verkko-thic/",
        verkko_hic_uri="s3://human-pangenomics/submissions/0624338D-1A6F-4E29-A276-D2C247FE0558--verkko-v2.1_intermediate_asms/{sm}/verkko-hi-c/",
    conda:
        "../envs/nucflag_bmk.yaml"
    shell:
        """
        aws s3 --no-sign-request sync {params.verkko_thic_uri} {output} \
        --include "*.bz*" \
        --exclude "*analysis/*" --exclude "*ribotin/*" || true

        # Separate naming convention.
        aws s3 --no-sign-request sync {params.verkko_hic_uri} {output} \
        --include "*.bz*" \
        --exclude "*analysis/*" --exclude "*ribotin/*" || true
        """


# Call to regenerate.
rule format_sample_sheet:
    input:
        script="workflow/rules/format_sample_list.py",
        sample_sheet="config/b2.tsv",
    output:
        formatted_sample_sheet=READS_SAMPLE_SHEET,
    conda:
        "../envs/nucflag_bmk.yaml"
    shell:
        """
        python {input.script} {input.sample_sheet} > {output}
        """


rule download_hifi:
    output:
        directory(os.path.join(READS_INPUT_DIR, "{sm}")),
    threads: 20
    params:
        files=lambda wc: SAMPLE_READ_URIS[str(wc.sm)],
    conda:
        "../envs/nucflag_bmk.yaml"
    shell:
        """
        mkdir -p {output}
        parallel -j {threads} aws s3 --no-sign-request cp {{}} {output}" ::: {params.files}
        """


rule generate_hifi_fofn:
    input:
        hifi_dir=rules.download_hifi.output,
    output:
        hifi_fofn=os.path.join(READS_INPUT_DIR, "{sm}.fofn"),
    shell:
        """
        find {input.hifi_dir} -type f > {output.hifi_fofn}
        """


checkpoint download_data_all:
    input:
        expand(rules.download_assemblies.output, sm=SAMPLE_NAMES),
        expand(rules.generate_hifi_fofn.output, sm=SAMPLE_NAMES),
    output:
        touch("download_data.done"),


rule download_bmk_data_all:
    input:
        rules.download_data_all.output,
    default_target: True
