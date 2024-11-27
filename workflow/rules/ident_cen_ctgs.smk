# Identify centromeric regions
# Extract centromeric regions


include: "utils.smk"
include: "common.smk"


# Convert rustybam stats bedfile by adjusting start and end positions.
rule format_hor_ref_aln_cen_contigs:
    input:
        aln_bed=os.path.join("results", f"{REF_NAME}_cens", "bed", "{sm}.bed"),
    output:
        cen_regions=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            "interm",
            "{sm}_cens.bed",
        ),
    conda:
        "../envs/tools.yaml"
    log:
        "logs/ident_cen_ctgs/format_hor_ref_aln_cen_contigs_{sm}.log",
    shell:
        # 1. chr2:91797392-95576642
        # 2. 3054999
        # 3. 3779251
        # 4. 3779251
        # 5. -
        # 6. h1tg000001l#1-110442987
        # 7. 15032098
        # 8. 15756783
        # 9. 110442987
        # 10. 99.97791
        # 11. 99.97032
        # 12. 99.89915
        # 13. 724023
        # 14. 160
        # 15. 27
        # 16. 28
        # 17. 69
        # 18. 502
        """
        awk -v OFS="\\t" 'NR>1 {{
            # Find starts/ends in contig name.
            match($1, ":(.+)-", starts);
            # Remove coords from ctg name
            gsub(":.*-.*", "", $1)
            # Print columns.
            print $1, $2+starts[1], $3+starts[1], $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18
        }}' {input} > {output} 2> {log}
        """


# Add unique regions in monomeric pq arms and intersect with alignments.
rule intersect_with_pq_arm:
    input:
        aln_cens_bed=rules.format_hor_ref_aln_cen_contigs.output,
        ref_monomeric_bed=config["ident_cen_ctgs"]["ref_cens_monomeric_regions"],
    output:
        qarms_cen_regions=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            "interm",
            "{sm}_pqarm_cens.bed",
        ),
    conda:
        "../envs/tools.yaml"
    log:
        "logs/ident_cen_ctgs/intersect_ref_cen_pqarm_{sm}.log",
    shell:
        """
        bedtools intersect -loj -a {input.aln_cens_bed} -b  {input.ref_monomeric_bed} > {output.qarms_cen_regions} 2> {log}
        """


# Map each centromeric contig to a chromosome based on intersection with CHM13 monomeric regions.
# Also attempt to reorient.
# No length threshold at this point.
rule map_collapse_cens:
    input:
        script="workflow/scripts/map_cens.py",
        regions=rules.intersect_with_pq_arm.output,
    output:
        cens_key=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            "interm",
            "{sm}_mapped_cens.bed",
        ),
        # old_name, new_name, coords, sample, chrom, is_reverse
        renamed_cens_key=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            "interm",
            "{sm}_renamed_cens.tsv",
        ),
    params:
        thr_ctg_len=0,
    conda:
        "../envs/py.yaml"
    log:
        "logs/ident_cen_ctgs/map_collapse_cens_{sm}.log",
    shell:
        """
        python {input.script} -i {input.regions} -t {params.thr_ctg_len} > {output.cens_key} 2> {log}
        awk -v OFS="\\t" '{{
            st=$2+1
            end=$3
            coords=st"-"end
            old_name=$1
            new_name="{wildcards.sm}_"$4"_"$1
            print old_name,new_name,coords,"{wildcards.sm}",$4,$6
        }}' {output.cens_key} > {output.renamed_cens_key} 2> {log}
        """


# Extract centromeric contigs.
use rule extract_and_index_fa as extract_cens_regions with:
    input:
        bed=rules.map_collapse_cens.output,
        fa=os.path.join(config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fa"),
    output:
        seq=temp(
            os.path.join(
                config["ident_cen_ctgs"]["output_dir"],
                "seq",
                "interm",
                "{sm}_centromeric_regions.fa",
            )
        ),
        idx=temp(
            os.path.join(
                config["ident_cen_ctgs"]["output_dir"],
                "seq",
                "interm",
                "{sm}_centromeric_regions.fa.fai",
            )
        ),
    log:
        "logs/ident_cen_ctgs/extract_regions_{sm}.log",
    conda:
        "../envs/tools.yaml"


rule ident_cen_ctgs_all:
    input:
        expand(
            rules.map_collapse_cens.output,
            sm=SAMPLE_NAMES,
        ),
