# Extract HOR arrays from reference.
rule extract_ref_hor_arrays:
    input:
        ref=config["align_asm_to_ref"]["config"]["ref"][REF_NAME],
        cens_regions=config["ident_cen_ctgs"]["ref_cens_500kbp_regions"],
    output:
        seq=os.path.join(
            config["extract_ref_hor_arrays"]["output_dir"],
            f"{REF_NAME}.hor_arrays.500kbp.fa",
        ),
        idx=os.path.join(
            config["extract_ref_hor_arrays"]["output_dir"],
            f"{REF_NAME}.hor_arrays.500kbp.fa.fai",
        ),
    log:
        "logs/extract_ref_hor_arrays.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        seqtk subseq {input.ref} {input.cens_regions} > {output.seq} 2> {log}
        samtools faidx {output.seq} &> {log}
        """


rule extract_ref_hor_arrays_all:
    input:
        rules.extract_ref_hor_arrays.output,
