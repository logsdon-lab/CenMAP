include: "common.smk"
include: "1-concat_asm.smk"
include: "4-ident_cen_ctgs.smk"


DNA_BRNN_OUTDIR = join(OUTPUT_DIR, "5-dna_brnn")
DNA_BRNN_LOGDIR = join(LOG_DIR, "5-dna_brnn")
DNA_BRNN_BMKDIR = join(BMK_DIR, "5-dna_brnn")


# Then split them into unique files per contig.
checkpoint split_fa_dnabrnn:
    input:
        fa=rules.extract_cens_regions.output.seq,
        rename_key=rules.map_collapse_cens.output.renamed_cens_key,
    output:
        directory(join(DNA_BRNN_OUTDIR, "seq", "interm", "{sm}")),
    log:
        join(DNA_BRNN_LOGDIR, "split_fa_{sm}.log"),
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
        join(DNA_BRNN_OUTDIR, "dna-nn", "dna-brnn"),
    conda:
        "../envs/dna-brnn.yaml"
    log:
        join(DNA_BRNN_LOGDIR, "compile_dna_brnn.log"),
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
        seqs=join(
            rules.split_fa_dnabrnn.output[0],
            "{fname}.fa",
        ),
    output:
        repeat_regions=join(
            DNA_BRNN_OUTDIR,
            "bed",
            "{sm}",
            "{fname}.bed",
        ),
    params:
        bin_dnabrnn=lambda wc, input: (
            input.bin_dnabrnn if input.bin_dnabrnn else "dna-brnn"
        ),
    threads: config["dna_brnn"]["threads"]
    resources:
        mem=config["dna_brnn"].get("mem", "4GB"),
    log:
        join(DNA_BRNN_LOGDIR, "dna_brnn_{sm}_{fname}.log"),
    benchmark:
        join(DNA_BRNN_BMKDIR, "dna_brnn_{sm}_{fname}.tsv")
    singularity:
        "docker://logsdonlab/dna-nn:latest"
    shell:
        """
        {params.bin_dnabrnn} -t {threads} -Ai {input.model} {input.seqs} > {output} 2> {log}
        """


# Filter dna-brnn output to alpha-satellite repeats only.
# Also trim to major asat HOR arrays.
rule filter_dnabrnn_sample_cens_regions:
    input:
        script=workflow.source_path("../scripts/filter_dnabrnn_output.py"),
        repeats=rules.run_dna_brnn.output,
        thresholds=config["dna_brnn"]["thr_file"],
    output:
        tmp_alr_ctgs=temp(
            join(
                DNA_BRNN_OUTDIR,
                "bed",
                "{sm}_filtered",
                "{fname}.bed",
            )
        ),
    params:
        chrom_name=lambda wc: get_chrom_name(wc.fname),
        repeat_type_filter=2,
    log:
        join(DNA_BRNN_LOGDIR, "filter_dnabrnn_{sm}_{fname}_cens_regions.log"),
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} \
        -i {input.repeats} \
        -o {output} \
        -t {input.thresholds} \
        --chr {params.chrom_name} \
        --repeat_type {params.repeat_type_filter} 2> {log}
        """


# TODO: Trim N's that are not directly adjacent to alpha-sat and update coordinates.
# seqkit locate -i results/dna_brnn/seq/interm/HG00171/HG00171_chrX_haplotype2-0000138\:54987333-61097711.fa -r -G -M -p "N+" | tail -n+2 | cut -f 1,5,6 | uniq
"""
# Find N's
seqkit locate -i {input.fa} -r -G -M -p "N+" | tail -n+2 | cut -f 1,5,6 | uniq > N_coords.bed
# relative coords to absolute coords
awk '{split($1, ctg_coords, ":"); split(ctg_coords[2], coords, "-"); print $1, coords[1] + $2, coords[1] + $3 }
# Check no intersect.
bedtools intersect
# Then if empty, adjust region start
"""


def dna_brnn_output(wc):
    outdir = checkpoints.split_fa_dnabrnn.get(sm=wc.sm).output[0]
    wcs = glob_wildcards(join(outdir, "{fname}.fa"))

    return expand(
        rules.filter_dnabrnn_sample_cens_regions.output,
        sm=wc.sm,
        fname=wcs.fname,
    )


rule aggregate_dnabrnn_alr_regions_by_sample:
    input:
        sample_cens=dna_brnn_output,
    output:
        join(
            DNA_BRNN_OUTDIR,
            "bed",
            "{sm}_all.bed",
        ),
    log:
        join(DNA_BRNN_LOGDIR, "aggregate_dnabrnn_alr_regions_by_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        cat {input.sample_cens} > {output} 2> {log}
        """


rule extract_alr_regions_by_sample:
    input:
        fa=rules.concat_asm.output.fa,
        bed=rules.aggregate_dnabrnn_alr_regions_by_sample.output,
    output:
        seq=temp(
            join(
                DNA_BRNN_OUTDIR,
                "seq",
                "interm",
                "{sm}_contigs.ALR.fa",
            )
        ),
        idx=temp(
            join(
                DNA_BRNN_OUTDIR,
                "seq",
                "interm",
                "{sm}_contigs.ALR.fa.fai",
            )
        ),
    params:
        **params_shell_extract_and_index_fa,
    log:
        join(DNA_BRNN_LOGDIR, "extract_alr_region_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        shell_extract_and_index_fa


rule dna_brnn_all:
    input:
        expand(
            rules.aggregate_dnabrnn_alr_regions_by_sample.output,
            sm=SAMPLE_NAMES,
        ),
    default_target: True
