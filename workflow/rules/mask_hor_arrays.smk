# Mask HOR arrays in reference genome T2T-CHM13
rule mask_hor_arrays:
    input:
        ref=config["align_asm_to_ref"]["config"]["ref"][REF_NAME],
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


# Extract masked HOR arrays from reference.
rule extract_masked_hor_arrays:
    input:
        masked_ref=rules.mask_hor_arrays.output,
        cens_regions=config["align_asm_to_ref"]["cens_500kbp_regions"],
    output:
        seq=os.path.join(
            config["mask_hor_arrays"]["output_dir"],
            f"{REF_NAME}.hor_arrays_masked.500kbp.fa",
        ),
        idx=os.path.join(
            config["mask_hor_arrays"]["output_dir"],
            f"{REF_NAME}.hor_arrays_masked.500kbp.fa.fai",
        ),
    log:
        "logs/extract_masked_hor_arrays.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        seqtk subseq {input.masked_ref} {input.cens_regions} > {output.seq} 2> {log}
        samtools faidx {output.seq} &> {log}
        """


# Extract HOR arrays from reference.
rule extract_hor_arrays:
    input:
        ref=config["align_asm_to_ref"]["config"]["ref"][REF_NAME],
        cens_regions=config["align_asm_to_ref"]["cens_500kbp_regions"],
    output:
        seq=os.path.join(
            config["mask_hor_arrays"]["output_dir"],
            f"{REF_NAME}.hor_arrays.500kbp.fa",
        ),
        idx=os.path.join(
            config["mask_hor_arrays"]["output_dir"],
            f"{REF_NAME}.hor_arrays.500kbp.fa.fai",
        ),
    log:
        "logs/extract_hor_arrays.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        seqtk subseq {input.ref} {input.cens_regions} > {output.seq} 2> {log}
        samtools faidx {output.seq} &> {log}
        """


rule extract_hor_arrays_all:
    input:
        rules.mask_hor_arrays.output,
        rules.extract_masked_hor_arrays.output,
        rules.extract_hor_arrays.output,
