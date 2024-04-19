REF_CENS_EDGE_LEN = round((500_000 + config["extract_ref_hor_arrays"].get("added_bases", 0)) / 1000)

rule adjust_ref_hor_arrays:
    input:
        config["ident_cen_ctgs"]["ref_cens_500kbp_regions"],
    output:
        os.path.join(
            config["extract_ref_hor_arrays"]["output_dir"],
            f"{REF_NAME}.hor_arrays.{REF_CENS_EDGE_LEN}kbp.bed",
        ),
    params:
        added_bases=config["extract_ref_hor_arrays"].get("added_bases", 0),
    shell:
        """
        awk -v OFS="\\t" '{{ print $1, $2-{params.added_bases}, $3+{params.added_bases} }}' {input} > {output}
        """


# Extract HOR arrays from reference.
rule extract_ref_hor_arrays:
    input:
        ref=REF_FA,
        cens_regions=rules.adjust_ref_hor_arrays.output,
    output:
        seq=os.path.join(
            config["extract_ref_hor_arrays"]["output_dir"],
            f"{REF_NAME}.hor_arrays.{REF_CENS_EDGE_LEN}kbp.fa",
        ),
        idx=os.path.join(
            config["extract_ref_hor_arrays"]["output_dir"],
            f"{REF_NAME}.hor_arrays.{REF_CENS_EDGE_LEN}kbp.fa.fai",
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
