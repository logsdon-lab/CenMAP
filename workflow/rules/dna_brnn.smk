import os


# https://github.com/lh3/dna-nn/tree/master
rule run_dna_brnn:
    input:
        model=config["dna_brnn"]["model"],
        seqs=rules.cens_rename_oriented_ctgs.output,
    output:
        repeat_regions=os.path.join(
            config["dna_brnn"]["output_dir"],
            "{sm}_centromeric_regions.renamed.{ort}.bed",
        ),
    threads: 20
    log:
        "logs/dna_brnn_{ort}_{sm}.log",
    benchmark:
        "benchmarks/dna_brnn_{ort}_{sm}.tsv"
    # No conda recipe. Use Dockerfile if not installed locally.
    singularity:
        "docker://koisland/hgsvc3:latest"
    shell:
        """
        dna-brnn -t {threads} -Ai {input.model} {input.seqs} > {output} 2> {log}
        """


# Run dna-brnn only on the ref centromeres as baseline expectation of ALR repeats.
# (Per chr)
# TODO: Replace with static values?
# /net/eichler/vol27/projects/hgsvc/nobackups/analysis/glennis/map_align_t2t_chm13/results/t2t_chm13_v2.0_cens/bed/centromeres/dna-brnn/chm13_cens.trimmed.bed
# dna-brnn on f"chm13.hor_arrays_masked.500kbp.fa"?
use rule run_dna_brnn as run_dna_brnn_ref_cens with:
    input:
        model=config["dna_brnn"]["model"],
        seqs=rules.extract_masked_hor_arrays.output,
    output:
        cens=os.path.join(
            config["dna_brnn"]["output_dir"], f"{REF_NAME}_cens.trimmed.bed"
        ),
    log:
        f"logs/dna_brnn_{REF_NAME}_cens.log",
    benchmark:
        f"benchmarks/dna_brnn_{REF_NAME}_cens.tsv"


# TODO: Script different from notebook?
# grep "chr1:" chm13_cens.trimmed.bed | \
# sed 's/:/\t/g' | sed 's/-/\t/g' | \
# awk -v OFS="\t" '{print $1, $2+$4, $2+$5, $6, $5-$4}' | awk '$4==2' | awk '$5>1000' > chr1_tmp.fwd.bed
rule filter_dnabrnn_ref_cens_regions:
    input:
        repeats=rules.run_dna_brnn_ref_cens.output,
    output:
        temp(os.path.join(config["dna_brnn"]["output_dir"], "{chr}_tmp.fwd.bed")),
    params:
        repeat_type_filter=2,
        repeat_len_thr=1000,
    log:
        "logs/filter_dnabrnn_ref_{chr}_cens_regions.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ grep "{wildcards.chr}:" {input.repeats} | \
        awk -v OFS="\\t" '{{print $1, $2, $3, $4, $3-$2}}' | \
        awk '$4=={params.repeat_type_filter} && $5>{params.repeat_len_thr}';}} > {output} 2> {log}
        """


# ex. >HG00171_chr16_haplotype1-0000003:4-8430174
# (Per chr + sample)
# grep "chr2_" ${sample}.renamed.fwd.bed | \
# sed 's/:/\t/g' | \
# sed 's/-/\t/g' | \
# awk -v OFS="\t" '{print $1"-"$2, $3+$5, $3+$6, $7, $6-$5}' | \
# awk '$4==2' | \
# awk '$5>1000' >> chr2_tmp.fwd.bed
rule filter_dnabrnn_sample_cens_regions:
    input:
        script="workflow/scripts/filter_cen_ctgs.py",
        repeats=rules.run_dna_brnn.output,
    output:
        temp(
            os.path.join(
                config["dna_brnn"]["output_dir"],
                "{chr}_{sm}_contigs.{ort}.ALR.bed",
            )
        ),
    params:
        split_cols=" ".join(["ctg_label", "ctg_num", "ctg_start", "ctg_stop"]),
        repeat_type_filter=2,
        repeat_len_thr=1000,
        is_forward_ort=lambda wc: "--forward" if wc.ort == "fwd" else "",
    log:
        "logs/filter_dnabrnn_{ort}_{sm}_{chr}_cens_regions.log",
    conda:
        "../env/py.yaml"
    shell:
        """
        python {input.script} filtdnabrnn \
        -i {input.repeats} \
        -o {output} \
        -c {wildcards.chr} \
        {params.is_forward_ort} \
        --columns_split {params.split_cols} \
        --repeat_type {params.repeat_type_filter} \
        --repeat_gt_length {params.repeat_len_thr} &> {log}
        """


# Calculate the start and end *-terms
# ex. chr2
# awk -v OFS="\t" '{print $1, $2-(*467987), $3+(*522450), $3-$2}' | \
rule get_dnabrnn_ref_cens_pos:
    input:
        script="workflow/scripts/get_cen_pos.py",
        filt_repeats=rules.filter_dnabrnn_ref_cens_regions.output,
    output:
        os.path.join(config["dna_brnn"]["output_dir"], "{chr}_cens_data.json"),
    log:
        "logs/get_dnabrnn_ref_{chr}_cens_data.log",
    conda:
        "../env/py.yaml"
    shell:
        """
        python {input.script} -i {input.filt_repeats} -o {output} &> {log}
        """


# TODO: Annotate
# /net/eichler/vol28/home/glogsdon/utilities/bedminmax.py \
# -i chr2_tmp.fwd.bed | \
# awk -v OFS="\t" '{print $1, $2, $3, $3-$2}' | \
# awk -v OFS="\t" '{print $1, $2-467987, $3+522450, $3-$2}' | \
# awk '$4>1000000' | \
# awk -v OFS="\t" '$2<0 {$2=0}1' > chr2_contigs.fwd.repeat.bed
rule aggregate_dnabrnn_alr_regions_by_chr:
    input:
        script="workflow/scripts/filter_cen_ctgs.py",
        cen_pos=rules.get_dnabrnn_ref_cens_pos.output,
        added_ref_cens=lambda wc: (
            rules.filter_dnabrnn_ref_cens_regions.output if wc.ort == "fwd" else []
        ),
        sample_cens=expand(
            os.path.join(
                config["dna_brnn"]["output_dir"],
                "{{chr}}_{sm}_contigs.{{ort}}.ALR.bed",
            ),
            sm=SAMPLE_NAMES,
        ),
    output:
        os.path.join(
            config["dna_brnn"]["output_dir"],
            "{chr}_contigs.{ort}.ALR.bed",
        ),
    params:
        repeat_len_thr=1_000_000,
        io_cols=" ".join(
            [
                "chr",
                "start",
                "end",
                "repeat_type",
                "repeat_length",
            ]
        ),
        grp_cols=" ".join(["chr", "repeat_type", "repeat_length"]),
        sort_cols=" ".join(["chr", "start"]),
    log:
        "logs/aggregate_dnabrnn_alr_regions_by_{chr}_{ort}.log",
    conda:
        "../env/py.yaml"
    shell:
        """
        start=$(jq .start {input.cen_pos})
        end=$(jq .end {input.cen_pos})

        # Aggregate and bedminmax all.
        # Select cols and calculate length.
        # Calculate vals.
        # Take only repeats greater than some value.
        # Take abs value.
        {{ python {input.script} bedminmax \
            -i {input.sample_cens} {input.added_ref_cens} \
            -ci {params.io_cols} \
            -co {params.io_cols} \
            -g {params.grp_cols} \
            -s {params.sort_cols} | \
        awk -v OFS="\\t" '{{print $1, $2, $3, $3-$2}}' | \
        awk -v START=$start END=$end OFS="\\t" '{{print $1, $2-start, $3+end, $3-$2}}' | \
        awk '$4>{params.repeat_len_thr}' | \
        awk -v OFS="\\t" '$2<0 {{$2=0}}1';}} > {output} 2> {log}
        """


rule dna_brnn_all:
    input:
        expand(rules.run_dna_brnn.output, sm=SAMPLE_NAMES, ort=ORIENTATION),
        rules.run_dna_brnn_ref_cens.output,
        expand(
            rules.filter_dnabrnn_ref_cens_regions.output,
            chr=CHROMOSOMES,
        ),
        expand(
            rules.filter_dnabrnn_sample_cens_regions.output,
            sm=SAMPLE_NAMES,
            ort=ORIENTATION,
            chr=CHROMOSOMES,
        ),
        expand(
            rules.get_dnabrnn_ref_cens_pos.output,
            chr=CHROMOSOMES,
        ),
        expand(
            rules.aggregate_dnabrnn_alr_regions_by_chr.output,
            ort=ORIENTATION,
            chr=CHROMOSOMES,
        ),
