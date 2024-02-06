# Align assemblies to reference.
# Pull alignment workflow from github.
# See output files here. Most of everything gets dumped in subdirs from results/{ref}/{ftype}/{asm}.{ftype}:
# https://github.com/mrvollger/asm-to-reference-alignment/blob/main/workflow/rules/reference_alignment.smk
module align_asm_to_ref:
    snakefile:
        github(
            "koisland/asm-to-reference-alignment",
            # "mrvollger/asm-to-reference-alignment",
            path="workflow/Snakefile",
            branch="remove_sm_num_index",
        )
    config:
        config["align_asm_to_ref"]["config"]


use rule * from align_asm_to_ref exclude dipcall, ideogram as asm_ref_*


# Override rules from ^ because input getter functions run before concat_asm so cannot detect file.
use rule alignment from align_asm_to_ref as asm_ref_alignment with:
    input:
        ref=lambda wc: config["align_asm_to_ref"]["config"]["ref"][wc.ref],
        # Take concat asm with all types included.
        query=rules.concat_asm.output,


use rule alignment2 from align_asm_to_ref as asm_ref_alignment2 with:
    input:
        ref_fasta=lambda wc: config["align_asm_to_ref"]["config"]["ref"][wc.ref],
        query=rules.concat_asm.output,
        # Weird. Blank if passing reference rule output from above. Use string instead.
        aln="temp/{ref}/{sm}.bam",


use rule all from align_asm_to_ref as align_asm_ref_all with:
    input:
        ra=rules.asm_ref_reference_alignment.input,
        SafFire=expand(
            rules.asm_ref_SafFire.output,
            sm=SAMPLE_NAMES,
            ref=config["align_asm_to_ref"]["config"]["ref"].keys(),
        ),
