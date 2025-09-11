include: "common.smk"
include: "2-concat_asm.smk"


SRF_OUTDIR = join(OUTPUT_DIR, "3-srf")
SRF_LOGDIR = join(LOG_DIR, "3-srf")
SRF_BMKDIR = join(BMK_DIR, "3-srf")
MONOMER_PERIODS = sorted(config["ident_cen_ctgs"]["monomer_periods"])


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
                    "input_file": rules.concat_asm.output.fa,
                    "parameters": {
                        "kmer_size": (
                            MONOMER_PERIODS[0] + 1
                            if MONOMER_PERIODS[0] % 2 == 0
                            else MONOMER_PERIODS[0]
                        ),
                        "exclude_kmers_lt_n": config["ident_cen_ctgs"][
                            "exclude_kmers_lt_n"
                        ],
                        "mon_periods": [*MONOMER_PERIODS, 42],
                        "perc_mon_len_diff": config["ident_cen_ctgs"][
                            "perc_mon_len_diff"
                        ]
                        / 100.0,
                    },
                }
                for sm in SAMPLE_NAMES
            },
            "workflow_dir": "workflow/rules/Snakemake-srf/workflow",
        }


use rule * from srf_sm as srf_*


checkpoint extract_filter_monomers:
    input:
        paf=ancient(rules.srf_map_motifs.output),
        monomers=ancient(rules.srf_get_monomers.output),
    output:
        bed_mon=join(
            SRF_OUTDIR,
            "bed",
            "{sm}_monomers_filtered.bed",
        ),
    params:
        # Look for alpha-sat and hsat-1a
        mon_periods=" ".join(str(e) for e in [*MONOMER_PERIODS, 42]),
        perc_mon_len_diff=config["ident_cen_ctgs"]["perc_mon_len_diff"] / 100.0,
    log:
        join(SRF_LOGDIR, "create_region_bed_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        srf-n-trf monomers \
        -p <(zcat {input.paf}) \
        -m {input.monomers} \
        -s {params.mon_periods} \
        -d {params.perc_mon_len_diff} > {output.bed_mon} 2> {log}
        """


rule merge_slop_region_bed:
    input:
        bed=rules.extract_filter_monomers.output.bed_mon,
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
        bp_merge=config["ident_cen_ctgs"]["bp_merge"],
        bp_min_length=config["ident_cen_ctgs"]["bp_min_length"],
        # But only use regions with alpha-sat.
        req_mon_periods=" ".join(str(e) for e in MONOMER_PERIODS),
        perc_mon_len_diff=config["ident_cen_ctgs"]["perc_mon_len_diff"] / 100.0,
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ srf-n-trf regions \
            -b {input.bed} \
            -m {params.bp_min_length} \
            -d {params.bp_merge} \
            -s {params.req_mon_periods} \
            --diff {params.perc_mon_len_diff} | \
        bedtools slop -i - -g {input.idx} -b {params.bp_slop} ;}} > {output} 2> {log}
        """


rule srf_all:
    input:
        expand(
            rules.merge_slop_region_bed.output,
            sm=SAMPLE_NAMES,
        ),
    default_target: True
