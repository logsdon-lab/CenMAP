include: "common.smk"
include: "2-concat_asm.smk"


SRF_OUTDIR = join(OUTPUT_DIR, "3-srf")
SRF_LOGDIR = join(LOG_DIR, "3-srf")
SRF_BMKDIR = join(BMK_DIR, "3-srf")


# Then split them into unique files per contig.
rule split_fa_srf:
    input:
        fa=rules.concat_asm.output.fa,
    output:
        temp(directory(join(SRF_OUTDIR, "seq", "interm", "{sm}"))),
    log:
        join(SRF_LOGDIR, "split_fa_{sm}.log"),
    params:
        split_dir=lambda wc, output: output[0],
    shell:
        """
        mkdir -p {params.split_dir}
        awk -v PIPE="|" '{{
            if (substr($0, 1, 1)==">") {{
                fname=substr($0,2)
                filename=("{params.split_dir}/" fname ".fa")
            }}
            # Ignore filenames with pipes.
            if (!index(filename, PIPE)) {{
                print $0 > filename
            }}
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
        bed=ancient(rules.srf_merge_files.output.bed),
        monomers=ancient(rules.srf_merge_files.output.monomers),
    output:
        bed=pipe(
            join(
                SRF_OUTDIR,
                "bed",
                "{sm}_interm.bed",
            )
        ),
    params:
        bp_group=config["ident_cen_ctgs"]["bp_group"],
        bp_min_length=config["ident_cen_ctgs"]["bp_min_length"],
        bp_merge=config["ident_cen_ctgs"]["bp_merge"],
        perc_mon_len_diff=config["ident_cen_ctgs"]["perc_mon_len_diff"],
    log:
        join(SRF_LOGDIR, "create_region_bed_{sm}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} \
        -i {input.bed} \
        -m {input.monomers} \
        -g {params.bp_group} \
        -l {params.bp_min_length} \
        -d {params.perc_mon_len_diff} \
        --merge_by {params.bp_merge} > {output.bed} 2> {log}
        """


rule slop_region_bed:
    input:
        bed=rules.create_region_bed.output,
        idx=rules.concat_asm.output.idx,
    output:
        join(
            SRF_OUTDIR,
            "bed",
            "{sm}.bed",
        ),
    log:
        join(SRF_LOGDIR, "format_region_bed_{sm}.log"),
    params:
        bp_slop=config["ident_cen_ctgs"]["bp_slop"],
    conda:
        "../envs/tools.yaml"
    shell:
        """
        bedtools slop -i {input.bed} -g {input.idx} -b {params.bp_slop} > {output} 2> {log}
        """


rule srf_all:
    input:
        expand(
            rules.slop_region_bed.output,
            sm=SAMPLE_NAMES,
        ),
    default_target: True
