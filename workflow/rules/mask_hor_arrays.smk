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
        rules.extract_hor_arrays.output,
