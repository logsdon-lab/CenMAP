
# Identify centromeric regions
# Extract centromeric regions

# TODO: Replace with wildcard for assembly
rule intersect_cent_regions:
    input:
        left="GM18989.bed",
        right="cenSat_Annotations_HORs.maxmin.v2.0.500kbp.bed"
    output:
        "GM18989_cens.bed"
    params:
        # Extra params
        extra=""
    log:
        "logs/intersect_cent.log"
    wrapper:
        "v3.3.3/bio/bedtools/intersect"

# TODO: This can be parallelized by adding unique fnames with wildcards.
# TODO: What is diff between this and next step?
# TODO: (low priority) bedmin_max.py uses pandas so we might be able to remove awk calls.
# * I'd have to write some tests to check behaviour is correct.
rule collapse_cens_contigs:
    input:
        cens_regions=rules.intersect_cent_regions.output
    output:
        "GM18989_1_centromeric_contigs.bed"
    shell:
        """
        awk -v OFS="\t" '{print $6, $7, $8, $1, $5}' {input.cens_regions} > tmp1.bed 

        ./bedminmax.py -i tmp1.bed -o tmp2.bed 

        awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $3-$2}' tmp2.bed | sort -k1,1 > {output} 
        """

# Reuse rule for vvv
# On output of asm_to_ref_alignment
use rule collapse_cens_contigs as collapse_cens_contigs_2 with:
    input:
        ""
    output:
        ""

# Reorient the regions.
rule reorient_regions:
    input:

    output:

    shell:
    

# Run DNA-BRNN.