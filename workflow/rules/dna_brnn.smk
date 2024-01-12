# https://github.com/lh3/dna-nn/tree/master
rule dna_brnn_fwd:
    input:
        model=config["dna_brnn"]["model"],
        cens=rules.rename_cens_fwd_ctgs.output,
    output:
        alr_regions="{}_centromeric_regions.renamed.fwd.bed",
    threads: 20
    # Requires installing dna-brnn. No conda recipe.
    envmodules:
        "dna-nn/0.1",
    shell:
        """
        dna-brnn -t {threads} -Ai {input.model} {input.cens} > {output}
        """


# TODO: Glennis has script to filter valid predict ALR per chr on UW cluster.
