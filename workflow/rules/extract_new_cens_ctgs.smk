# TODO: This might be redone to go at start. Wait until later.


rule extract_new_oriented_cens_regions:
    input:
        regions=rules.run_dna_brnn.output.repeat_regions,
        combined_assembly=os.path.join(
            config["ident_cen_ctgs"]["comb_assemblies_dir"],
            "{sm}.vrk-ps-sseq.asm-comb-dedup.fa.gz",
        ),
    output:
        os.path.join(config["dna_brnn"]["output_dir"], "{sm}_{num}_contigs.{ort}.fa"),
    wildcard_constraints:
        ort="fwd|rev",
    params:
        added_cmds=lambda wc: "" if wc.ort == "fwd" else "| seqtk seq -r",
    log:
        "logs/extract_new_{ort}_cens_regions_{sm}_{num}.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        seqtk subseq {input.combined_assembly} {input.regions} {params.added_cmds} > {output}
        """


# HG00171_chr4_haplotype1-0000002:1892469-12648706        293621  293950  1
# |HG00171|chr4|haplotype1-0000002|1892469-12648706|293621|293950|1
RENAME_NEW_CTGS_CFG = {
    # TODO: Is this correct? This is before processed.
    "bed_input_regions": rules.run_dna_brnn.output.repeat_regions,
    "fa_assembly": rules.extract_new_oriented_cens_regions.output,
    "output_dir": os.path.join(config["dna_brnn"]["output_dir"], "new_cens"),
    "samples": SAMPLE_NAMES,
    "log_dir": "logs/rename_cens",
    "bed_find_col": 3,
    "bed_replace_w_joined_cols": (1, 2, 3),
}


module rename_new_cens_ctgs:
    snakefile:
        "rename_ctgs.smk"
    config:
        RENAME_NEW_CTGS_CFG


use rule * from rename_new_cens_ctgs as new_cens_*
