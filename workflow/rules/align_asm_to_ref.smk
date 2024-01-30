# Add additional ref to be aligned to all samples.
config["align_asm_to_ref"]["config"]["ref"][
    f"{REF_NAME}_cens_hardmasked"
] = os.path.join(
    config["mask_hor_arrays"]["output_dir"], f"{REF_NAME}.hor_arrays_masked.fa"
)

config["align_asm_to_ref"]["config"]["ref"][f"{REF_NAME}_cens"] = os.path.join(
    config["mask_hor_arrays"]["output_dir"],
    f"{REF_NAME}.hor_arrays_masked.500kbp.fa",
)


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


use rule * from align_asm_to_ref as asm_ref_*


# Override rules from ^ because input getter functions run before concat_asm so cannot detect file.
use rule alignment from align_asm_to_ref as asm_ref_alignment with:
    input:
        ref=config["align_asm_to_ref"]["config"]["ref"][REF_NAME],
        # Take concat asm with all types included.
        query=rules.concat_asm.output,


use rule alignment2 from align_asm_to_ref as asm_ref_alignment2 with:
    input:
        ref_fasta=config["align_asm_to_ref"]["config"]["ref"][REF_NAME],
        query=rules.concat_asm.output,
        # Weird. Blank if passing reference rule output from above. Use string instead.
        aln="temp/{ref}/{sm}.bam",
