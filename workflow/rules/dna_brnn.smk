import os


# https://github.com/lh3/dna-nn/tree/master
rule run_dna_brnn:
    input:
        model=config["dna_brnn"]["model"],
        seqs=rules.rename_cens_oriented_ctgs.output,
    output:
        repeat_regions=os.path.join(
            config["dna_brnn"]["output_dir"],
            "{sm}_centromeric_regions.renamed.{ort}.bed",
        ),
    threads: 20
    log:
        "logs/dna_brnn_{ort}_{sm}.log",
    # No conda recipe. Use Dockerfile if not installed locally.
    singularity:
        "docker://koisland/hgsvc3:latest"
    shell:
        """
        dna-brnn -t {threads} -Ai {input.model} {input.seqs} > {output} 2> {log}
        """


# Run dna-brnn only on the ref centromeres as baseline expectation of ALR repeats.
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


# grep "chr1:" chm13_cens.trimmed.bed | \
# sed 's/:/\t/g' | sed 's/-/\t/g' | \
# awk -v OFS="\t" '{print $1, $2+$4, $2+$5, $6, $5-$4}' | awk '$4==2' | awk '$5>1000' > chr1_tmp.fwd.bed
rule filter_ref_cens_regions:
    input:
        script="workflow/scripts/filter_cen_ctgs.py",
        cens=rules.run_dna_brnn_ref_cens.output,
    output:
        temp(os.path.join(config["dna_brnn"]["output_dir"], "{chr}_tmp.fwd.bed")),
    params:
        # Columns created when splitting name.
        split_cols=" ".join(["ctg_label", "ctg_start", "ctg_stop"]),
        repeat_type_filter=2,
        repeat_len_thr=1000,
    conda:
        "../env/py.yaml"
    shell:
        """
        python {input.script} filtdnabrnn \
        -i {input.cens} \
        -o {output} \
        -c {wildcards.chr} \
        --columns_split {params.split_cols} \
        --repeat_type {params.repeat_type_filter} \
        --repeat_len_thr {params.repeat_len_thr}
        """


# TODO: Glennis has script to filter valid predict ALR per chr on UW cluster.
rule filter_sample_cens_regions:
    input:
        lambda wc: rules.run_dna_brnn.output,
    output:
        temp("{chr}_{sm}_contigs.{ort}.ALR.bed"),
    shell:
        """
        """


def get_chr_stats(wc):
    pass


rule aggregate_alr_by_chr:
    input:
        lambda wc: expand(
            rules.filter_sample_cens_regions,
            chr=[wc.chr],
            ort=[wc.ort],
            sm=SAMPLES_DF.index,
        ),
    output:
        "{chr}_contigs.{ort}.ALR.bed",
    shell:
        """
        """


rule dna_brnn_all:
    input:
        expand(rules.run_dna_brnn.output, sm=SAMPLES_DF.index, ort=ORIENTATION),
        expand(
            rules.filter_sample_cens_regions.output,
            sm=SAMPLES_DF.index,
            ort=ORIENTATION,
            chr=CHROMOSOMES,
        ),
        rules.run_dna_brnn_ref_cens.output,
        # ...
