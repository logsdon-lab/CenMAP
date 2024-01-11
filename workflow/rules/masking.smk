# Hardmask centromeres


# Mask HOR arrays in reference genome T2T-CHM13
rule mask_hor_arrays:
    input:
        ref="/net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/T2T-CHM13v2.fasta",
        hor_arrays="/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/t2t_chr8/ape_cens/chm1/asat_annotation/chm13_v2.0_hor_arrays.bed",
    output:
        masked_ref="T2T-CHM13v2.hor_arrays_masked.fa",
    conda:
        "env/masking.yaml"
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
        rules.mask_hor_arrays.output,
        "cenSat_Annotations_HORs.maxmin.v2.0.500kbp.bed",
    output:
        "T2T-CHM13v2.hor_arrays_masked.500kbp.fa",
    log:
        "logs/extract_masked_hor_arrays.log",
    params:
        command="subseq",
        extra="",
    wrapper:
        "v3.3.3/bio/seqtk"

# Reindex masked reference genome.
rule index_masked_ref:
    input:
        rules.extract_masked_hor_arrays.output
    output:
        "T2T-CHM13v2.hor_arrays_masked.500kbp.fai",
    log:
        "logs/index_masked_ref.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v3.3.3/bio/samtools/index"