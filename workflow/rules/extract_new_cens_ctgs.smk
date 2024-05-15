# HG00171_chr22_h2tg000021l#1-94189493:1993043-3660421
# HG00171_chr22_h2tg000021l#1-94189493
rule create_format_orient_cens_list:
    input:
        regions=rules.run_dna_brnn.output.repeat_regions,
    output:
        os.path.join(
            config["new_cens"]["output_dir"],
            "bed",
            "{sm}_contigs.{ort}.list",
        ),
    log:
        "logs/extract_new_cens_ctgs/format_{ort}_cens_list_{sm}.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{  sed -e 's/:/\\t/g' {input} | cut -f 1 | sort | uniq;}} > {output} 2> {log}
        """


use rule extract_and_index_fa as extract_new_oriented_cens_regions with:
    input:
        bed=rules.create_format_orient_cens_list.output,
        fa=rules.asm_rename_ctgs.output,
    output:
        seq=os.path.join(
            config["new_cens"]["output_dir"], "seq", "{sm}_contigs.{ort}.fa"
        ),
        idx=os.path.join(
            config["new_cens"]["output_dir"], "seq", "{sm}_contigs.{ort}.fa.fai"
        ),
    params:
        added_cmds=lambda wc: "" if wc.ort == "fwd" else "| seqtk seq -r",
    log:
        "logs/extract_new_cens_ctgs/extract_new_{ort}_cens_regions_{sm}.log",


# Corresponds to
# /net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/hgsvc3/...
# ...map_align_t2t_chm13/results/t2t_chm13_v2.0_cens/bed/centromeres/dna-brnn/centromeric_regions/extract_ALR_regions.bash
use rule extract_and_index_fa as extract_alr_region_ref_by_chr with:
    input:
        fa=REF_FA,
        bed=lambda wc: expand(
            rules.aggregate_dnabrnn_alr_regions_by_chr.output,
            chr=[wc.chr],
            ort=["fwd"],
        ),
    output:
        seq=temp(
            os.path.join(
                config["new_cens"]["output_dir"],
                "seq",
                f"{REF_NAME}_{{chr}}_contigs.fwd.ALR.fa",
            )
        ),
        idx=temp(
            os.path.join(
                config["new_cens"]["output_dir"],
                "seq",
                f"{REF_NAME}_{{chr}}_contigs.fwd.ALR.fa.fai",
            )
        ),
    params:
        added_cmds="",
    log:
        f"logs/extract_new_cens_ctgs/extract_alr_region_{REF_NAME}_{{chr}}.log",


use rule extract_and_index_fa as extract_alr_region_sample_by_chr with:
    input:
        fa=rules.extract_new_oriented_cens_regions.output.seq,
        bed=rules.aggregate_dnabrnn_alr_regions_by_chr.output,
    output:
        seq=temp(
            os.path.join(
                config["new_cens"]["output_dir"],
                "seq",
                "{chr}_{sm}_contigs.{ort}.ALR.fa",
            )
        ),
        idx=temp(
            os.path.join(
                config["new_cens"]["output_dir"],
                "seq",
                "{chr}_{sm}_contigs.{ort}.ALR.fa.fai",
            )
        ),
    log:
        "logs/extract_new_cens_ctgs/extract_alr_region_{sm}_{chr}_{ort}.log",


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
        seq=os.path.join(
            config["new_cens"]["output_dir"], "seq", "{chr}_contigs.{ort}.ALR.fa"
        ),
        idx=os.path.join(
            config["new_cens"]["output_dir"], "seq", "{chr}_contigs.{ort}.ALR.fa.fai"
        ),
    conda:
        "../env/tools.yaml"
    log:
        "logs/extract_new_cens_ctgs/merge_alr_regions_by_{chr}_{ort}.log",
    shell:
        """
        cat {input.ref_ctgs_fa} {input.sm_ctgs_fa} > {output.seq} 2> {log}
        cat {input.ref_ctgs_fai} {input.sm_ctgs_fai} > {output.idx} 2> {log}
        """


# ex. NA12329 chrX    haplotype1      0000017 92036218        98212716        6176499 4617952 6176499 6176500
# ex. HG00171 chr22   h1tg000027l     1       26260313        1       4719177 4719177 48      4719177 4719178
rule make_new_cens_bed_file:
    input:
        script="workflow/scripts/filter_cen_ctgs.py",
        faidx=expand(
            rules.merge_alr_regions_by_chr.output.idx,
            ort=ORIENTATION,
            sm=SAMPLE_NAMES,
            chr=CHROMOSOMES,
        ),
    output:
        tmp_fmt_alr_bed=temp(
            os.path.join(
                config["new_cens"]["output_dir"],
                "bed",
                "fmt_{sm}_ALR_regions.bed",
            )
        ),
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
        uniq;}} > {output.tmp_fmt_alr_bed} 2> {log}

        {{ python {input.script} bedminmax \
            -i {output.tmp_fmt_alr_bed} \
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
            ort=ORIENTATION,
            chr=CHROMOSOMES,
        ),
        expand(rules.make_new_cens_bed_file.output, sm=SAMPLE_NAMES),
