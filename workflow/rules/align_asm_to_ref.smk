ALN_CFG = {
    "ref": {f"{REF_NAME}_cens": rules.extract_ref_hor_arrays.output.seq},
    "sm": {
        sm: os.path.join(config["concat_asm"]["output_dir"], f"{sm}-asm-comb-dedup.fa")
        for sm in SAMPLE_NAMES
    },
    "temp_dir": os.path.join(config["align_asm_to_ref"]["output_dir"], "temp"),
    "output_dir": config["align_asm_to_ref"]["output_dir"],
    "logs_dir": "logs/align_asm_to_ref",
    "benchmarks_dir": "benchmarks/align_asm_to_ref",
    "aln_threads": config["align_asm_to_ref"]["threads"],
    "mm2_opts": "-x asm20 --secondary=no -s 25000 -K 15G",
}


# Align assemblies to reference.
# Pull alignment workflow from github.
# See output files here. Most of everything gets dumped in subdirs from results/{ref}/{ftype}/{asm}.{ftype}:
# https://github.com/mrvollger/asm-to-reference-alignment/blob/main/workflow/rules/reference_alignment.smk
module align_asm_to_ref:
    snakefile:
        github(
            "koisland/asm-to-reference-alignment",
            path="workflow/Snakefile",
            branch="minimal",
        )
    config:
        ALN_CFG


use rule * from align_asm_to_ref as asm_ref_*


use rule all from align_asm_to_ref as align_asm_ref_all with:
    input:
        rules.asm_ref_all.input,
