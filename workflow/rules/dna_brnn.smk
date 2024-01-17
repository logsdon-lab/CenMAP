# https://github.com/lh3/dna-nn/tree/master
rule dna_brnn_fwd:
    input:
        model=config["dna_brnn"]["model"],
        cens=rules.rename_cens_fwd_ctgs.output,
    output:
        alr_regions="{sm}_centromeric_regions.renamed.fwd.bed",
    threads: 20
    log:
        "logs/dna_brnn_fwd_{sm}.log",
    # No conda recipe. Use Dockerfile if not installed locally.
    singularity:
        "docker://koisland/hgsvc3:latest"
    shell:
        """
        dna-brnn -t {threads} -Ai {input.model} {input.cens} > {output} 2> {log}
        """


use rule dna_brnn_fwd as dna_brnn_rev with:
    input:
        model=config["dna_brnn"]["model"],
        cens=rules.rename_cens_rev_ctgs.output,
    output:
        alr_regions="{sm}_centromeric_regions.renamed.rev.bed",
    log:
        "logs/dna_brnn_rev_{sm}.log",


# TODO: Glennis has script to filter valid predict ALR per chr on UW cluster.


rule dna_brnn_all:
    input:
        expand(rules.dna_brnn_fwd.output, sm=SAMPLES_DF.index),
        expand(rules.dna_brnn_rev.output, sm=SAMPLES_DF.index),
        # ...
