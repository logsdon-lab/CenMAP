# Hardmask centromeres


# Mask HOR arrays in reference genome T2T-CHM13
rule mask_hor_arrays:
    input:
        ref="/net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/T2T-CHM13v2.fasta",
        hor_arrays="/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/t2t_chr8/ape_cens/chm1/asat_annotation/chm13_v2.0_hor_arrays.bed",
    output:
        masked_ref="T2T-CHM13v2.hor_arrays_masked.fa",
    conda:
        "env/env.yaml"
    log:
        "logs/mask_hor_array_ref.log",
    shell:
        """
        bedtools maskfasta \
        -fi {input.ref} \
        -bed {input.hor_arrays} \
        -fo {output.masked_ref} &> {log}
        """


# Extract masked HOR arrays from reference?
rule extract_masked_hor_arrays:
    input:
        masked_ref=rules.mask_hor_arrays.output,
        hor_array_regions="cenSat_Annotations_HORs.maxmin.v2.0.500kbp.bed",
    output:
        "T2T-CHM13v2.hor_arrays_masked.500kbp.fa",
    log:
        "logs/extract_masked_hor_arrays.log",
    conda:
        "env/env.yaml"
    shell:
        """
        seqtk subseq {input.masked_ref} {input.hor_array_regions} > {output} 2> {log}
        """

# Then index them.
rule index_masked_hor_array:
    input:
        rules.extract_masked_hor_arrays.output,
    output:
        "T2T-CHM13v2.hor_arrays_masked.500kbp.fai",
    log:
        "logs/index_masked_hor_array.log",
    shell:
        """
        samtools faidx {input} 2> {log}
        """
