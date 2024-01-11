
# Hardmask centromeres


rule mask_hor_arrays:
    input:
        ref = "/net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/T2T-CHM13v2.fasta",
        bed = "/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/t2t_chr8/ape_cens/chm1/asat_annotation/chm13_v2.0_hor_arrays.bed"
    output:
        masked_ref = "T2T-CHM13v2.hor_arrays_masked.fa"
    shell:
        """
        bedtools maskfasta \
        -fi {input.ref} \
        -bed {input.bed} \
        -fo {output.masked_ref}
        """

# rule extract_masked_hor_arrays:
#     input: