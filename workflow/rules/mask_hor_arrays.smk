import os


REF_MASK_OUTPUT = []

# Mask HOR arrays in reference genome T2T-CHM13
if "cens_hardmasked" in config["mask_hor_arrays"]["added_alignments"]:

    rule mask_hor_arrays:
        input:
            ref=config["align_asm_to_ref"]["config"]["ref"][REF_NAME],
            # TODO: Are the annotations per chr or all grouped together?
            asat_annotations=config["align_asm_to_ref"]["asat_annotations"],
        output:
            os.path.join(
                config["mask_hor_arrays"]["output_dir"],
                f"{REF_NAME}.hor_arrays_masked.fa",
            ),
        conda:
            "../env/tools.yaml"
        log:
            "logs/mask_hor_array_ref.log",
        shell:
            """
            bedtools maskfasta \
            -fi {input.ref} \
            -bed {input.asat_annotations} \
            -fo {output} &> {log}
            """

    REF_MASK_OUTPUT.append(rules.mask_hor_arrays.output)


if "cens" in config["mask_hor_arrays"]["added_alignments"]:

    # Extract masked HOR arrays from reference?
    rule extract_masked_hor_arrays:
        input:
            masked_ref=rules.mask_hor_arrays.output,
            cens_regions=config["align_asm_to_ref"]["cens_500kbp_regions"],
        output:
            os.path.join(
                config["mask_hor_arrays"]["output_dir"],
                f"{REF_NAME}.hor_arrays_masked.500kbp.fa",
            ),
        log:
            "logs/extract_masked_hor_arrays.log",
        conda:
            "../env/tools.yaml"
        shell:
            """
            seqtk subseq {input.masked_ref} {input.cens_regions} > {output} 2> {log}
            """

    # Then index them.
    rule index_masked_hor_array:
        input:
            rules.extract_masked_hor_arrays.output,
        output:
            os.path.join(
                config["mask_hor_arrays"]["output_dir"],
                f"{REF_NAME}.hor_arrays_masked.500kbp.fa.fai",
            ),
        conda:
            "../env/tools.yaml"
        log:
            "logs/index_masked_hor_array.log",
        shell:
            """
            samtools faidx {input} &> {log}
            """

    REF_MASK_OUTPUT.append(rules.extract_masked_hor_arrays.output)
    REF_MASK_OUTPUT.append(rules.index_masked_hor_array.output)


rule mask_hor_arrays_all:
    input:
        REF_MASK_OUTPUT,
