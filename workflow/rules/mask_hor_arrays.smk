# Hardmask centromeres


# Mask HOR arrays in reference genome T2T-CHM13
rule mask_hor_arrays:
    input:
        ref=config["asm_to_ref"]["config"]["ref"]["CHM13"],
        # TODO: Are the annotations per chr or all grouped together?
        asat_annotations=config["asm_to_ref"]["asat_annotations"]
    output:
        masked_ref="data/T2T-CHM13v2.hor_arrays_masked.fa",
    conda:
        "../env/tools.yaml"
    log:
        "logs/mask_hor_array_ref.log",
    shell:
        """
        bedtools maskfasta \
        -fi {input.ref} \
        -bed {input.asat_annotations} \
        -fo {output.masked_ref} &> {log}
        """


# Extract masked HOR arrays from reference?
rule extract_masked_hor_arrays:
    input:
        masked_ref=rules.mask_hor_arrays.output,
        cens_regions=config["asm_to_ref"]["cens_500kbp_regions"],
    output:
        "data/T2T-CHM13v2.hor_arrays_masked.500kbp.fa",
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
        "data/T2T-CHM13v2.hor_arrays_masked.500kbp.fai",
    conda:
        "../env/tools.yaml"
    log:
        "logs/index_masked_hor_array.log",
    shell:
        """
        samtools faidx {input} 2> {log}
        """
