# Pull alignment workflow from github
# See output files here. Most of everything gets dumped in subdirs from results/{ref}/{ftype}/{asm}.{ftype}:
# https://github.com/mrvollger/asm-to-reference-alignment/blob/main/workflow/rules/reference_alignment.smk
module asm_to_ref_align:
    snakefile:
        github(
            "mrvollger/asm-to-reference-alignment",
            path="workflow/Snakefile",
            tag="v0.1",
        )
    config: config

use rule * from asm_to_ref_align as ref_align_*
