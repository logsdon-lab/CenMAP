include: "common.smk"


rule compile_dna_brnn:
    output:
        os.path.join(config["dna_brnn"]["output_dir"], "dna-nn", "dna-brnn"),
    conda:
        "../envs/dna-brnn.yaml"
    log:
        "logs/dna_brnn/compile_dna_brnn.log",
    params:
        url="https://github.com/lh3/dna-nn",
        tmp_log="compile_dna_brnn.log",
        output_dir=lambda wc, output: os.path.dirname(str(output)),
    shell:
        """
        log_path=$(realpath {log})
        rm -rf {params.output_dir} && git clone {params.url} {params.output_dir} 2> "${{log_path}}"
        cd {params.output_dir}
        conda_env=$(which gcc | sed 's/\\/bin\\/gcc//g')
        C_INCLUDE_PATH="${{conda_env}}/include/" make >> "${{log_path}}" 2>> "${{log_path}}"
        """


# https://github.com/lh3/dna-nn/tree/master
rule run_dna_brnn:
    input:
        bin_dnabrnn=rules.compile_dna_brnn.output if not IS_CONTAINERIZE_CMD else [],
        model=config["dna_brnn"]["model"],
        seqs=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "seq",
            "{sm}_centromeric_regions.fa",
        ),
    output:
        repeat_regions=os.path.join(
            config["dna_brnn"]["output_dir"],
            "{sm}",
            "{sm}_centromeric_regions.renamed.bed",
        ),
    params:
        bin_dnabrnn=lambda wc, input: (
            input.bin_dnabrnn if input.bin_dnabrnn else "dna-brnn"
        ),
    threads: config["dna_brnn"]["threads"]
    resources:
        mem=config["dna_brnn"].get("mem", "4GB"),
    log:
        "logs/dna_brnn/dna_brnn_{sm}.log",
    benchmark:
        "benchmarks/dna_brnn/dna_brnn_{sm}.tsv"
    singularity:
        "docker://logsdonlab/dna-nn:latest"
    shell:
        """
        {params.bin_dnabrnn} -t {threads} -Ai {input.model} {input.seqs} > {output} 2> {log}
        """


# Run dna-brnn only on the ref centromeres as baseline expectation of ALR repeats.
# (Per chr)
# /net/eichler/vol27/projects/hgsvc/nobackups/analysis/glennis/map_align_t2t_chm13/results/t2t_chm13_v2.0_cens/bed/centromeres/dna-brnn/chm13_cens.trimmed.bed
# dna-brnn on f"chm13.hor_arrays_masked.500kbp.fa"?
use rule run_dna_brnn as run_dna_brnn_ref_cens with:
    input:
        model=config["dna_brnn"]["model"],
        seqs=os.path.join(
            config["extract_ref_hor_arrays"]["output_dir"],
            f"{REF_NAME}.hor_arrays.{REF_CENS_EDGE_LEN}kbp.bed",
        ),
    output:
        cens=os.path.join(
            config["dna_brnn"]["output_dir"], REF_NAME, f"{REF_NAME}_cens.trimmed.bed"
        ),
    log:
        f"logs/dna_brnn/dna_brnn_{REF_NAME}_cens.log",
    benchmark:
        f"benchmarks/dna_brnn/dna_brnn_{REF_NAME}_cens.tsv"


# grep "chr1:" chm13_cens.trimmed.bed | \
# sed 's/:/\t/g' | sed 's/-/\t/g' | \
# awk -v OFS="\t" '{print $1, $2+$4, $2+$5, $6, $5-$4}' | awk '$4==2' | awk '$5>1000' > chr1_tmp.fwd.bed
rule filter_dnabrnn_ref_cens_regions:
    input:
        script="workflow/scripts/filter_dnabrnn_output.py",
        repeats=(
            rules.run_dna_brnn_ref_cens.output
            if config["dna_brnn"].get("ref_alr_file") is None
            else config["dna_brnn"]["ref_alr_file"]
        ),
        thresholds=config["dna_brnn"]["thr_file"],
    output:
        temp(
            os.path.join(
                config["dna_brnn"]["output_dir"], REF_NAME, "{chr}_tmp.fwd.bed"
            )
        ),
    params:
        repeat_type_filter=2,
    log:
        "logs/dna_brnn/filter_dnabrnn_ref_{chr}_cens_regions.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} \
        -i {input.repeats} \
        -o {output} \
        -t {input.thresholds} \
        --chr {wildcards.chr} \
        --repeat_type {params.repeat_type_filter} 2> {log}
        """


# ex. >HG00171_chr19_h1tg000004l#1-58112442
# ex. >HG00171_chr16_haplotype1-0000003:4-8430174
# HG00171_chr16_haplotype1|0000003|4|8430174|start|end|type
# (Per chr + sample)
# grep "chr2_" ${sample}.renamed.fwd.bed | \
# sed 's/:/\t/g' | \
# sed 's/-/\t/g' | \
# awk -v OFS="\t" '{print $1"-"$2, $3+$5, $3+$6, $7, $6-$5}' | \
# awk '$4==2' | \
# awk '$5>1000' >> chr2_tmp.fwd.bed
use rule filter_dnabrnn_ref_cens_regions as filter_dnabrnn_sample_cens_regions with:
    input:
        script="workflow/scripts/filter_dnabrnn_output.py",
        repeats=rules.run_dna_brnn.output,
        thresholds=config["dna_brnn"]["thr_file"],
    output:
        tmp_alr_ctgs=temp(
            os.path.join(
                config["dna_brnn"]["output_dir"],
                "{sm}",
                "{chr}_{sm}_contigs.ALR.bed",
            )
        ),
    params:
        repeat_type_filter=2,
    log:
        "logs/dna_brnn/filter_dnabrnn_{sm}_{chr}_cens_regions.log",
    conda:
        "../envs/py.yaml"


# /net/eichler/vol28/home/glogsdon/utilities/bedminmax.py (modified bedminmax) \
# -i chr2_tmp.fwd.bed | \
# awk -v OFS="\t" '{print $1, $2, $3, $3-$2}' | \
# awk -v OFS="\t" '{print $1, $2-467987, $3+522450, $3-$2}' | \
# awk '$4>1000000' | \
# awk -v OFS="\t" '$2<0 {$2=0}1' > chr2_contigs.fwd.repeat.bed
rule aggregate_dnabrnn_alr_regions_by_chr:
    input:
        added_ref_cens=rules.filter_dnabrnn_ref_cens_regions.output,
        sample_cens=lambda wc: expand(
            rules.filter_dnabrnn_sample_cens_regions.output,
            sm=SAMPLE_NAMES,
            chr=[wc.chr],
        ),
    output:
        os.path.join(
            config["dna_brnn"]["output_dir"],
            "{chr}_contigs.ALR.bed",
        ),
    params:
        repeat_len_thr=dnabrnn_alr_region_threshold,
        # "chr", "start", "end", "repeat_type", "repeat_len"
        grp_cols="1,4",
        # After grouping
        # "chr", "repeat_type", "start", "end"
        sort_cols="-k1 -k3n",
        bp_edges=500_000,
        # dna-brnn output may not contain repeats from a chr.
        allow_empty="--allow_empty",
    log:
        "logs/dna_brnn/aggregate_dnabrnn_alr_regions_by_{chr}.log",
    conda:
        "../envs/tools.yaml"
    # Aggregate and bedminmax all.
    # Select cols and calculate length.
    # Add 500 kbp on both ends.
    # Take only repeats greater than some value.
    # Take abs value.
    shell:
        """
        {{ bedtools groupby \
        -i <(cat {input.sample_cens} {input.added_ref_cens} | sort | uniq) \
        -g {params.grp_cols} \
        -c 2,3 \
        -o min,max | \
        sort {params.sort_cols} | \
        awk -v BP_EDGES={params.bp_edges} -v OFS="\\t" '{{
            len=$4-$3
            if (len > {params.repeat_len_thr}) {{
                new_start=$3-BP_EDGES
                new_end=$4+BP_EDGES
                new_start=((new_start < 0) ? 0 : new_start)
                print $1, new_start, new_end, $4-$3
            }}
        }}';}} > {output} 2> {log}
        """


rule dna_brnn_all:
    input:
        expand(
            rules.aggregate_dnabrnn_alr_regions_by_chr.output,
            chr=CHROMOSOMES,
        ),
