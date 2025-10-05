include: "utils.smk"
include: "common.smk"
include: "3-srf.smk"


ASM_REF_OUTDIR = join(OUTPUT_DIR, "4-align_asm_to_ref")
ASM_REF_LOGDIR = join(LOG_DIR, "4-align_asm_to_ref")
ASM_REF_BMKDIR = join(BMK_DIR, "4-align_asm_to_ref")

ALN_CFG = {
    "ref": (
        {
            REF_NAME: (
                config["ident_cen_ctgs"]["reference"]
                if config["ident_cen_ctgs"].get("reference")
                else REF_FA
            )
        }
        if CHROMOSOMES
        else {}
    ),
    "sm": {sm: rules.concat_asm.output.fa for sm in SAMPLE_NAMES},
    "temp_dir": join(ASM_REF_OUTDIR, "temp"),
    "output_dir": ASM_REF_OUTDIR,
    "logs_dir": ASM_REF_LOGDIR,
    "benchmarks_dir": ASM_REF_BMKDIR,
    "aln_threads": config["ident_cen_ctgs"]["threads_aln"],
    "aln_mem": config["ident_cen_ctgs"]["mem_aln"],
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
    default_target: True
