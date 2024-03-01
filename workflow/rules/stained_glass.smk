if config["stained_glass"]["input_dir"]:
    INPUT_FA_DIR = config["stained_glass"]["input_dir"]
else:
    INPUT_FA_DIR = config["humas_hmmer"]["input_dir"]


rule index_fa_for_stained_glass:
    input:
        fa=os.path.join(INPUT_FA_DIR, "{fname}.fa"),
    output:
        idx=os.path.join(INPUT_FA_DIR, "{fname}.fa.fai"),
    conda:
        "../env/stained_glass.yaml"
    log:
        "logs/index_fa_for_stained_glass+{fname}.log",
    shell:
        """
        # Index file.
        samtools faidx {input.fa} 2>> {log}
        """


rule run_stained_glass:
    input:
        snakefile="workflow/scripts/StainedGlass/workflow/Snakefile",
        fa=os.path.join(INPUT_FA_DIR, "{fname}.fa"),
        idx=rules.index_fa_for_stained_glass.output,
    output:
        directory(
            "{fname}"
            + f".{config['stained_glass']['window']}.{config['stained_glass']['mm_f']}_figures"
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
        # snakemake --configfile config/config.yaml --unlock 2> {log}

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
def stained_glass_outputs_no_input_dir(wc):
    _ = checkpoints.split_cens_for_humas_hmmer.get(**wc).output
    fnames = glob_wildcards(os.path.join(INPUT_FA_DIR, "{fname}.fa")).fname
    return expand(rules.run_stained_glass.output, fname=fnames)


# Conditionally change based on provided input dir.
if config["stained_glass"]["input_dir"] is None:

    rule stained_glass_all:
        input:
            stained_glass_outputs_no_input_dir,
        output:
            temp(touch("/tmp/stained_glass_{chr}.done")),

else:
    fnames = glob_wildcards(os.path.join(INPUT_FA_DIR, "{fname}.fa")).fname

    rule stained_glass_all:
        input:
            expand(rules.run_stained_glass.output, fname=fnames),
