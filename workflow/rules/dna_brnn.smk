import os


# https://github.com/lh3/dna-nn/tree/master
rule run_dna_brnn:
    input:
        model=config["dna_brnn"]["model"],
        cens=rules.rename_cens_oriented_ctgs.output,
    output:
        alr_regions=os.path.join(
            config["dna_brnn"]["output_dir"],
            "{sm}_centromeric_regions.renamed.{ort}.bed",
        ),
    threads: 20
    log:
        "logs/dna_brnn_{ort}_{sm}.log",
    # No conda recipe. Use Dockerfile if not installed locally.
    singularity:
        "docker://koisland/hgsvc3:latest"
    shell:
        """
        dna-brnn -t {threads} -Ai {input.model} {input.cens} > {output} 2> {log}
        """


# TODO: Glennis has script to filter valid predict ALR per chr on UW cluster.
rule filter_valid_alr:
    input:
        lambda wc: rules.run_dna_brnn.output,
    output:
        temp("{chr}_{sm}_contigs.{ort}.ALR.bed"),
    shell:
        """
        """


def get_chr_stats(wc):
    pass


rule aggregate_alr_by_chr:
    input:
        lambda wc: expand(
            rules.filter_valid_alr, chr=[wc.chr], ort=[wc.ort], sm=SAMPLES_DF.index
        ),
    output:
        "{chr}_contigs.{ort}.ALR.bed",
    shell:
        """
        """


rule dna_brnn_all:
    input:
        expand(rules.run_dna_brnn.output, sm=SAMPLES_DF.index, ort=ORIENTATION),
        expand(
            rules.filter_valid_alr.output,
            sm=SAMPLES_DF.index,
            ort=ORIENTATION,
            chr=CHROMOSOMES,
        ),
        # ...
