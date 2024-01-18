# Add additional ref to be aligned to all samples.
if "cens_hardmasked" in config["mask_hor_arrays"]["added_alignments"]:
    config["align_asm_to_ref"]["config"]["ref"][
        f"{REF_NAME}_cens_hardmasked"
    ] = os.path.join(
        config["mask_hor_arrays"]["output_dir"], f"{REF_NAME}.hor_arrays_masked.fa"
    )

if "cens" in config["mask_hor_arrays"]["added_alignments"]:
    config["align_asm_to_ref"]["config"]["ref"][f"{REF_NAME}_cens"] = os.path.join(
        config["mask_hor_arrays"]["output_dir"],
        f"{REF_NAME}.hor_arrays_masked.500kbp.fa",
    )


# Align assemblies to reference.
# Pull alignment workflow from github.
# * TODO: Still need to see original config.yaml.
# See output files here. Most of everything gets dumped in subdirs from results/{ref}/{ftype}/{asm}.{ftype}:
# https://github.com/mrvollger/asm-to-reference-alignment/blob/main/workflow/rules/reference_alignment.smk
module align_asm_to_ref:
    snakefile:
        github(
            "mrvollger/asm-to-reference-alignment",
            path="workflow/Snakefile",
            branch="main",
        )
    config:
        config["align_asm_to_ref"]["config"]


use rule * from align_asm_to_ref as ref_align_*
