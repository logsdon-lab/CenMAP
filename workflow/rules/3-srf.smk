include: "common.smk"
include: "1-concat_asm.smk"


SRF_OUTDIR = join(OUTPUT_DIR, "3-srf")
SRF_LOGDIR = join(LOG_DIR, "3-srf")
SRF_BMKDIR = join(BMK_DIR, "3-srf")


# Then split them into unique files per contig.
rule split_fa_srf:
    input:
        fa=rules.concat_asm.output.fa,
    output:
        directory(join(SRF_OUTDIR, "seq", "interm", "{sm}")),
    log:
        join(SRF_LOGDIR, "split_fa_{sm}.log"),
    params:
        split_dir=lambda wc, output: output[0],
    shell:
        """
        mkdir -p {params.split_dir}
        awk '{{
            if (substr($0, 1, 1)==">") {{
                fname=substr($0,2)
                filename=("{params.split_dir}/" fname ".fa")
            }}
            print $0 > filename
        }}' {input.fa} 2> {log}
        """


module srf_sm:
    snakefile:
        "Snakemake-srf/workflow/Snakefile"
    config:
        {
            "output_dir": SRF_OUTDIR,
            "log_dir": SRF_LOGDIR,
            "benchmark_dir": SRF_BMKDIR,
            "threads": config["ident_cen_ctgs"]["threads_srf"],
            "mem": config["ident_cen_ctgs"]["mem_srf"],
            "samples": {
                sm: {
                    "input_dir": expand(rules.split_fa_srf.output, sm=sm),
                    # Use default.
                    "parameters": {"kmer_size": 151, "exclude_kmers_lt_n": 3},
                }
                for sm in SAMPLE_NAMES
            },
        }


use rule * from srf_sm as srf_*


rule create_region_bed:
    input:
        script=workflow.source_path("../scripts/filter_srf_bed.py"),
        bed=ancient(rules.srf_merge_files_n_cleanup.output.bed),
        monomers=ancient(rules.srf_merge_files_n_cleanup.output.monomers),
    output:
        bed=join(
            SRF_OUTDIR,
            "bed",
            "{sm}.bed",
        ),
    log:
        join(SRF_LOGDIR, "format_region_bed_{sm}.log"),
    params:
        added_bases=config["ident_cen_ctgs"]["added_bases"],
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} -i {input.bed} -m {input.monomers} --bp_slop {params.added_bases} > {output.bed} 2> {log}
        """


rule extract_alr_regions_by_sample:
    input:
        fa=rules.concat_asm.output.fa,
        bed=rules.create_region_bed.output,
    output:
        seq=join(
            SRF_OUTDIR,
            "seq",
            "interm",
            "{sm}.fa",
        ),
        idx=join(
            SRF_OUTDIR,
            "seq",
            "interm",
            "{sm}.fa.fai",
        ),
    params:
        **params_shell_extract_and_index_fa,
    log:
        join(SRF_LOGDIR, "extract_alr_region_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        shell_extract_and_index_fa


rule srf_all:
    input:
        ancient(
            expand(
                rules.extract_alr_regions_by_sample.output,
                sm=SAMPLE_NAMES,
            )
        ),
    default_target: True
