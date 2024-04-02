NEW_CENS_OUTPUT_DIR = os.path.join(config["dna_brnn"]["output_dir"], "new_cens")


rule create_format_orient_cens_list:
    input:
        regions=rules.run_dna_brnn.output.repeat_regions,
    output:
        os.path.join(
            NEW_CENS_OUTPUT_DIR,
            "{sm}_contigs.{ort}.list",
        ),
    log:
        "logs/format_{ort}_cens_list_{sm}.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ cut -f 1 {input} | sort | uniq;}} > {output} 2> {log}
        """


use rule extract_and_index_fa as extract_new_oriented_cens_regions with:
    input:
        bed=rules.create_format_orient_cens_list.output,
        fa=rules.asm_rename_ctgs.output,
    output:
        seq=os.path.join(NEW_CENS_OUTPUT_DIR, "{sm}_contigs.{ort}.fa"),
        idx=os.path.join(NEW_CENS_OUTPUT_DIR, "{sm}_contigs.{ort}.fa.fai"),
    params:
        added_cmds=lambda wc: "" if wc.ort == "fwd" else "| seqtk seq -r",
    log:
        "logs/extract_new_{ort}_cens_regions_{sm}.log",


# Corresponds to
# /net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/hgsvc3/...
# ...map_align_t2t_chm13/results/t2t_chm13_v2.0_cens/bed/centromeres/dna-brnn/centromeric_regions/extract_ALR_regions.bash
use rule extract_and_index_fa as extract_alr_region_ref_by_chr with:
    input:
        fa=config["align_asm_to_ref"]["config"]["ref"][REF_NAME],
        bed=lambda wc: expand(
            rules.aggregate_dnabrnn_alr_regions_by_chr.output,
            chr=[wc.chr],
            ort=["fwd"],
        ),
    output:
        seq=temp(
            os.path.join(NEW_CENS_OUTPUT_DIR, f"{REF_NAME}_{{chr}}_contigs.fwd.ALR.fa")
        ),
        idx=temp(
            os.path.join(
                NEW_CENS_OUTPUT_DIR, f"{REF_NAME}_{{chr}}_contigs.fwd.ALR.fa.fai"
            )
        ),
    params:
        added_cmds="",
    log:
        f"logs/extract_alr_region_{REF_NAME}_{{chr}}.log",


use rule extract_and_index_fa as extract_alr_region_sample_by_chr with:
    input:
        fa=rules.extract_new_oriented_cens_regions.output.seq,
        bed=rules.aggregate_dnabrnn_alr_regions_by_chr.output,
    output:
        seq=temp(os.path.join(NEW_CENS_OUTPUT_DIR, "{chr}_{sm}_contigs.{ort}.ALR.fa")),
        idx=temp(
            os.path.join(NEW_CENS_OUTPUT_DIR, "{chr}_{sm}_contigs.{ort}.ALR.fa.fai")
        ),
    log:
        "logs/extract_alr_region_{sm}_{chr}_{ort}.log",


rule merge_alr_regions_by_chr:
    input:
        sm_ctgs_fa=lambda wc: expand(
            rules.extract_alr_region_sample_by_chr.output.seq,
            chr=[wc.chr],
            sm=SAMPLE_NAMES,
            ort=[wc.ort],
        ),
        ref_ctgs_fa=rules.extract_alr_region_ref_by_chr.output.seq,
        sm_ctgs_fai=lambda wc: expand(
            rules.extract_alr_region_sample_by_chr.output.idx,
            chr=[wc.chr],
            sm=SAMPLE_NAMES,
            ort=[wc.ort],
        ),
        ref_ctgs_fai=lambda wc: rules.extract_alr_region_ref_by_chr.output.idx,
    output:
        seq=os.path.join(NEW_CENS_OUTPUT_DIR, "{chr}_contigs.{ort}.ALR.fa"),
        idx=os.path.join(NEW_CENS_OUTPUT_DIR, "{chr}_contigs.{ort}.ALR.fa.fai"),
    conda:
        "../env/tools.yaml"
    log:
        "logs/merge_alr_regions_by_{chr}_{ort}.log",
    shell:
        """
        cat {input.ref_ctgs_fa} {input.sm_ctgs_fa} > {output.seq} 2> {log}
        cat {input.ref_ctgs_fai} {input.sm_ctgs_fai} > {output.idx} 2> {log}
        """


rule extract_new_cens_all:
    input:
        expand(
            rules.create_format_orient_cens_list.output,
            sm=SAMPLE_NAMES,
            ort=ORIENTATION,
            chr=CHROMOSOMES,
        ),
        expand(
            rules.extract_new_oriented_cens_regions.output,
            sm=SAMPLE_NAMES,
            ort=ORIENTATION,
            chr=CHROMOSOMES,
        ),
        expand(
            rules.extract_alr_region_ref_by_chr.output,
            sm=SAMPLE_NAMES,
            ort=ORIENTATION,
            chr=CHROMOSOMES,
        ),
        expand(
            rules.extract_alr_region_sample_by_chr.output,
            sm=SAMPLE_NAMES,
            ort=ORIENTATION,
            chr=CHROMOSOMES,
        ),
        expand(
            rules.merge_alr_regions_by_chr.output,
            sm=SAMPLE_NAMES,
            ort=ORIENTATION,
            chr=CHROMOSOMES,
        ),
