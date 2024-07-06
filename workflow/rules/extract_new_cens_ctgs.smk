# Corresponds to
# /net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/hgsvc3/...
# ...map_align_t2t_chm13/results/t2t_chm13_v2.0_cens/bed/centromeres/dna-brnn/centromeric_regions/extract_ALR_regions.bash
use rule extract_and_index_fa as extract_alr_region_ref_by_chr with:
    input:
        fa=REF_FA,
        bed=lambda wc: expand(
            rules.aggregate_dnabrnn_alr_regions_by_chr.output,
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
    params:
        added_cmds="",
    log:
        f"logs/extract_new_cens_ctgs/extract_alr_region_{REF_NAME}_{{chr}}.log",


use rule extract_and_index_fa_w_rc_bed as extract_alr_region_sample_by_chr with:
    input:
        fa=rules.asm_rename_ctgs.output,
        bed=rules.aggregate_dnabrnn_alr_regions_by_chr.output,
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


rule merge_alr_regions_by_chr:
    input:
        sm_ctgs_fa=lambda wc: expand(
            rules.extract_alr_region_sample_by_chr.output.seq,
            chr=[wc.chr],
            sm=SAMPLE_NAMES,
        ),
        ref_ctgs_fa=rules.extract_alr_region_ref_by_chr.output.seq,
        sm_ctgs_fai=lambda wc: expand(
            rules.extract_alr_region_sample_by_chr.output.idx,
            chr=[wc.chr],
            sm=SAMPLE_NAMES,
        ),
        ref_ctgs_fai=lambda wc: rules.extract_alr_region_ref_by_chr.output.idx,
    output:
        seq=os.path.join(
            config["new_cens"]["output_dir"], "seq", "{chr}_contigs.ALR.fa"
        ),
        idx=os.path.join(
            config["new_cens"]["output_dir"], "seq", "{chr}_contigs.ALR.fa.fai"
        ),
    conda:
        "../env/tools.yaml"
    log:
        "logs/extract_new_cens_ctgs/merge_alr_regions_by_{chr}.log",
    shell:
        """
        cat {input.ref_ctgs_fa} {input.sm_ctgs_fa} > {output.seq} 2> {log}
        cat {input.ref_ctgs_fai} {input.sm_ctgs_fai} > {output.idx} 2> {log}
        """


rule fmt_new_cens_bed_file:
    input:
        faidx=expand(
            rules.merge_alr_regions_by_chr.output.idx,
            sm=SAMPLE_NAMES,
            chr=CHROMOSOMES,
        ),
    output:
        temp(
            os.path.join(
                config["new_cens"]["output_dir"],
                "bed",
                "fmt_{sm}_ALR_regions.bed",
            )
        ),
    conda:
        "../env/tools.yaml"
    log:
        "logs/extract_new_cens_ctgs/fmt_{sm}_new_cens_bed_file.log",
    shell:
        # Only filter for sample to avoid malformed output ref cols in alr bed.
        """
        {{ cat {input.faidx} | \
        awk -v OFS="\\t" '{{
            if ($1 ~ "^{wildcards.sm}") {{
                match($1, "^(.+):", ctgs);
                match($1, ":(.+)-", starts);
                match($1, ".*-(.+)$", ends);
                print ctgs[1], starts[1], ends[1]
            }}
        }}' | \
        sort | \
        uniq;}} > {output} 2> {log}
        """


# ex. NA12329 chrX    haplotype1      0000017 92036218        98212716        6176499 4617952 6176499 6176500
# ex. HG00171 chr22   h1tg000027l     1       26260313        1       4719177 4719177 48      4719177 4719178
rule make_new_cens_bed_file:
    input:
        script="workflow/scripts/filter_cen_ctgs.py",
        bed=rules.fmt_new_cens_bed_file.output,
    output:
        alr_bed=os.path.join(
            config["new_cens"]["output_dir"], "bed", "{sm}_ALR_regions.bed"
        ),
    params:
        io_cols=" ".join(["ctg", "start", "end"]),
        grp_cols=" ".join(["ctg"]),
        sort_cols=" ".join(["ctg", "start"]),
    conda:
        "../env/py.yaml"
    log:
        "logs/extract_new_cens_ctgs/make_{sm}_bed_files_for_plot.log",
    shell:
        """
        {{ python {input.script} bedminmax \
            -i {input.bed} \
            -ci {params.io_cols} \
            -co {params.io_cols} \
            -g {params.grp_cols} \
            -s {params.sort_cols} | \
        awk -v OFS="\\t" '{{ print $1, $2, $3, $3-$2, $4}}';}} > {output.alr_bed} 2>> {log}
        """


rule extract_new_cens_all:
    input:
        expand(
            rules.merge_alr_regions_by_chr.output,
            sm=SAMPLE_NAMES,
            chr=CHROMOSOMES,
        ),
        expand(rules.make_new_cens_bed_file.output, sm=SAMPLE_NAMES),
