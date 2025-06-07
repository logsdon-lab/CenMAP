include: "utils.smk"
include: "common.smk"
include: "2-extract_ref_hor_arrays.smk"
include: "3-srf.smk"


ASM_REF_OUTDIR = join(OUTPUT_DIR, "4-align_asm_to_ref")
ASM_REF_LOGDIR = join(LOG_DIR, "4-align_asm_to_ref")
ASM_REF_BMKDIR = join(BMK_DIR, "4-align_asm_to_ref")

ALN_CFG = {
    "ref": {f"{REF_NAME}_cens": rules.extract_ref_hor_arrays.output.seq},
    "sm": {sm: ancient(expand(rules.extract_alr_regions_by_sample.output.seq, sm=sm)) for sm in SAMPLE_NAMES},
    "temp_dir": join(ASM_REF_OUTDIR, "temp"),
    "output_dir": ASM_REF_OUTDIR,
    "logs_dir": ASM_REF_LOGDIR,
    "benchmarks_dir": ASM_REF_BMKDIR,
    "aln_threads": config["align_asm_to_ref"]["threads"],
    "aln_mem": config["align_asm_to_ref"]["mem"],
    "mm2_opts": "-x asm20 --secondary=no -s 25000 -K 8G",
}


# Align assemblies to reference.
# Pull alignment workflow from github.
# See output files here. Most of everything gets dumped in subdirs from results/{ref}/{ftype}/{asm}.{ftype}:
# https://github.com/mrvollger/asm-to-reference-alignment/blob/main/workflow/rules/reference_alignment.smk
module align_asm_to_ref:
    snakefile:
        "asm-to-reference-alignment/workflow/Snakefile"
    config:
        ALN_CFG


use rule * from align_asm_to_ref as asm_ref_*


use rule all from align_asm_to_ref as align_asm_ref_all with:
    input:
        rules.asm_ref_all.input,
