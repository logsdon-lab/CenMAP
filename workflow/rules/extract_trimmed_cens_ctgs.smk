include: "common.smk"


# Corresponds to
# /net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/hgsvc3/...
# ...map_align_t2t_chm13/results/t2t_chm13_v2.0_cens/bed/centromeres/dna-brnn/centromeric_regions/extract_ALR_regions.bash
use rule extract_and_index_fa as extract_alr_region_ref_by_chr with:
    input:
        fa=(
            config["align_asm_to_ref"]["reference"]
            if config["align_asm_to_ref"].get("reference")
            else REF_FA
        ),
        bed=lambda wc: expand(
            os.path.join(
                config["dna_brnn"]["output_dir"],
                "{chr}_contigs.ALR.bed",
            ),
            chr=[wc.chr],
        ),
    output:
        seq=temp(
            os.path.join(
                config["new_cens"]["output_dir"],
                "seq",
                f"{REF_NAME}_{{chr}}_contigs.ALR.fa",
            )
        ),
        idx=temp(
            os.path.join(
                config["new_cens"]["output_dir"],
                "seq",
                f"{REF_NAME}_{{chr}}_contigs.ALR.fa.fai",
            )
        ),
    log:
        f"logs/extract_new_cens_ctgs/extract_alr_region_{REF_NAME}_{{chr}}.log",


use rule extract_and_index_fa as extract_alr_region_sample_by_chr with:
    input:
        fa=os.path.join(config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fa"),
        bed=os.path.join(
            config["dna_brnn"]["output_dir"],
            "{chr}_contigs.ALR.bed",
        ),
    output:
        seq=temp(
            os.path.join(
                config["new_cens"]["output_dir"],
                "seq",
                "{chr}_{sm}_contigs.ALR.fa",
            )
        ),
        idx=temp(
            os.path.join(
                config["new_cens"]["output_dir"],
                "seq",
                "{chr}_{sm}_contigs.ALR.fa.fai",
            )
        ),
    log:
        "logs/extract_new_cens_ctgs/extract_alr_region_{sm}_{chr}.log",


rule make_new_cens_bed_file:
    input:
        faidx=lambda wc: expand(
            rules.extract_alr_region_sample_by_chr.output.idx,
            chr=CHROMOSOMES,
            sm=[wc.sm],
        ),
    output:
        alr_bed=os.path.join(
            config["new_cens"]["output_dir"], "bed", "{sm}_ALR_regions.bed"
        ),
    conda:
        "../envs/tools.yaml"
    log:
        "logs/extract_new_cens_ctgs/fmt_{sm}_new_cens_bed_file.log",
    shell:
        # Only filter for sample to avoid malformed output ref cols in alr bed.
        """
        {{ cat {input.faidx} | \
        awk -v OFS="\\t" '{{
            match($1, "^(.+):", ctgs);
            match($1, ":(.+)-", starts);
            match($1, ".*-(.+)$", ends);
            print ctgs[1], starts[1], ends[1]
        }}' | \
        sort | \
        uniq;}} > {output} 2> {log}
        """


rule extract_trimmed_cens_all:
    input:
        expand(rules.make_new_cens_bed_file.output, sm=SAMPLE_NAMES),
