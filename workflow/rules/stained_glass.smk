
rule run_stained_glass:
    input:
        snakefile="workflow/scripts/StainedGlass/workflow/Snakefile",
        fa=os.path.join(config["humas_hmmer"]["input_dir"], "{fname}.fa"),
    output:
        directory(
            f"{{fname}}.{config['stained_glass']['window']}.{config['stained_glass']['mm_f']}_figures"
        ),
    conda:
        "../env/stained_glass.yaml"
    threads: config["stained_glass"]["threads"]
    params:
        window=config["stained_glass"]["window"],
        mm_f=config["stained_glass"]["mm_f"],
        nbatch=4,
        target_rule="make_figures",
    log:
        "logs/stained_glass_{fname}.log",
    benchmark:
        "benchmarks/stained_glass_{fname}.txt"
    shell:
        """
        # Remove lock.
        snakemake --configfile config/config.yaml --unlock 2> {log}

        # Index file.
        samtools faidx {input.fa} 2>> {log}

        # Run StainedGlass
        snakemake -s workflow/scripts/StainedGlass/workflow/Snakefile \
        --use-conda \
        --config \
        sample="{wildcards.fname}" \
        fasta={input.fa} \
        window={params.window} \
        mm_f={params.mm_f} \
        nbatch={params.nbatch} \
        --cores {threads} \
        -p \
        {params.target_rule} 2>> {log}
        """


# https://stackoverflow.com/a/63040288
def stained_glass_outputs(wc):
    _ = checkpoints.split_cens_for_humas_hmmer.get(**wc).output
    fnames = glob_wildcards(
        os.path.join(config["humas_hmmer"]["input_dir"], "{fname}.fa")
    ).fname

    return expand(rules.run_stained_glass.output, fname=fnames)
# rule all:
#     input:
#         expand(
#             rules.run_stained_glass.output, fname=["NA19240_rc_chr9_haplotype2-0000100:10693477-13672460"]
#         )
#     default_target: True
