include: "common.smk"


# Then split them into unique files per contig.
checkpoint split_fa_dnabrnn:
    input:
        fa=rules.extract_cens_regions.output.seq,
        rename_key=rules.map_collapse_cens.output.renamed_cens_key,
    output:
        directory(
            os.path.join(config["dna_brnn"]["output_dir"], "seq", "interm", "{sm}")
        ),
    log:
        "logs/dna_brnn/split_fa_{sm}.log",
    params:
        split_dir=lambda wc, output: output[0],
    shell:
        """
        mkdir -p {params.split_dir}
        awk '{{
            # Read key values in first file.
            if (FNR == NR) {{
                # Add coords to name.
                kv[$1":"$3]=$2":"$3;
                next;
            }}

            if (substr($0, 1, 1)==">") {{
                fname=substr($0,2)
                repl_fname=kv[fname]
                filename=("{params.split_dir}/" repl_fname ".fa")
            }}
            print $0 > filename
        }}' {input.rename_key} {input.fa} 2> {log}
        """


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
        # Only compile dna-brnn if conda-only.
        bin_dnabrnn=(
            rules.compile_dna_brnn.output if IS_CONDA and not IS_SINGULARITY else []
        ),
        model=config["dna_brnn"]["model"],
        seqs=os.path.join(
            rules.split_fa_dnabrnn.output[0],
            "{sm}_{chr}_{fname}.fa",
        ),
    output:
        repeat_regions=os.path.join(
            config["dna_brnn"]["output_dir"],
            "bed",
            "{chr}",
            "{sm}_{chr}_{fname}_centromeric_regions.bed",
        ),
    params:
        bin_dnabrnn=lambda wc, input: (
            input.bin_dnabrnn if input.bin_dnabrnn else "dna-brnn"
        ),
    threads: config["dna_brnn"]["threads"]
    resources:
        mem=config["dna_brnn"].get("mem", "4GB"),
    log:
        "logs/dna_brnn/dna_brnn_{sm}_{chr}_{fname}.log",
    benchmark:
        "benchmarks/dna_brnn/dna_brnn_{sm}_{chr}_{fname}.tsv"
    singularity:
        "docker://logsdonlab/dna-nn:latest"
    shell:
        """
        {params.bin_dnabrnn} -t {threads} -Ai {input.model} {input.seqs} > {output} 2> {log}
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
rule filter_dnabrnn_sample_cens_regions:
    input:
        script="workflow/scripts/filter_dnabrnn_output.py",
        repeats=rules.run_dna_brnn.output,
        thresholds=config["dna_brnn"]["thr_file"],
    output:
        tmp_alr_ctgs=temp(
            os.path.join(
                config["dna_brnn"]["output_dir"],
                "bed",
                "{chr}",
                "{sm}_{chr}_{fname}.ALR.bed",
            )
        ),
    params:
        repeat_type_filter=2,
    log:
        "logs/dna_brnn/filter_dnabrnn_{sm}_{chr}_{fname}_cens_regions.log",
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


def dna_brnn_output(wc):
    outdir = checkpoints.split_fa_dnabrnn.get(sm=wc.sm).output[0]
    sms_all, fnames_all, chrs_all = [], [], []
    wcs = glob_wildcards(os.path.join(outdir, "{sm}_{chr}_{fname}.fa"))
    for sm, chrom, fname in zip(wcs.sm, wcs.chr, wcs.fname):
        sms_all.append(sm)
        chrs_all.append(chrom)
        fnames_all.append(fname)

    return expand(
        rules.filter_dnabrnn_sample_cens_regions.output,
        zip,
        sm=sms_all,
        chr=chrs_all,
        fname=fnames_all,
    )


# /net/eichler/vol28/home/glogsdon/utilities/bedminmax.py (modified bedminmax) \
# -i chr2_tmp.fwd.bed | \
# awk -v OFS="\t" '{print $1, $2, $3, $3-$2}' | \
# awk -v OFS="\t" '{print $1, $2-467987, $3+522450, $3-$2}' | \
# awk '$4>1000000' | \
# awk -v OFS="\t" '$2<0 {$2=0}1' > chr2_contigs.fwd.repeat.bed
rule aggregate_dnabrnn_alr_regions_by_chr:
    input:
        sample_cens=dna_brnn_output,
    output:
        os.path.join(
            config["dna_brnn"]["output_dir"],
            "bed",
            "{sm}_contigs.ALR.bed",
        ),
    log:
        "logs/dna_brnn/aggregate_dnabrnn_alr_regions_by_{sm}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        cat {input.sample_cens} > {output} 2> {log}
        """


rule dna_brnn_all:
    input:
        expand(
            rules.aggregate_dnabrnn_alr_regions_by_chr.output,
            sm=SAMPLE_NAMES,
        ),
