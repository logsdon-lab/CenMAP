
include: "common.smk"
include: "utils.smk"
include: "1-concat_asm.smk"
include: "7-fix_cens_w_repeatmasker.smk"


HUMAS_SD_OUTDIR = join(OUTPUT_DIR, "8-humas_sd")
HUMAS_SD_LOGDIR = join(LOG_DIR, "8-humas_sd")
HUMAS_SD_BMKDIR = join(BMK_DIR, "8-humas_sd")


rule extract_cens_for_humas_sd:
    input:
        fa=rules.rename_reort_asm.output.fa,
        bed=rules.make_complete_cens_bed.output,
    output:
        seq=temp(
            join(
                HUMAS_SD_OUTDIR,
                "seq",
                "interm",
                "{sm}_cens.fa",
            )
        ),
        idx=temp(
            join(
                HUMAS_SD_OUTDIR,
                "seq",
                "interm",
                "{sm}_cens.fa.fai",
            )
        ),
    log:
        join(HUMAS_SD_LOGDIR, "extract_cens_for_humas_sd_{sm}.log"),
    params:
        # Seqtk outputs 1-based coords which causes name issues.
        bed=lambda wc, input: f"""<(awk -v OFS="\\t" '{{print $1, $2-1, $3}}' {input.bed})""",
        added_cmds="",
    shell:
        shell_extract_and_index_fa


checkpoint split_cens_for_humas_sd:
    input:
        fa=rules.extract_cens_for_humas_sd.output.seq,
    output:
        touch(
            join(
                HUMAS_SD_OUTDIR,
                "split_cens_for_humas_sd_{sm}.done",
            )
        ),
    log:
        join(HUMAS_SD_LOGDIR, "split_{sm}_cens_for_humas_sd.log"),
    params:
        split_dir=HUMAS_CENS_SPLIT_DIR,
    conda:
        "../envs/tools.yaml"
    shell:
        # https://gist.github.com/astatham/621901
        """
        mkdir -p {params.split_dir}
        awk '{{
            if (substr($0, 1, 1)==">") {{
                filename=("{params.split_dir}/" substr($0,2) ".fa")
            }}
            print $0 > filename
        }}' {input.fa} 2> {log}
        """


module HumAS_SD:
    snakefile:
        "Snakemake-HumAS-SD/workflow/Snakefile"
    config:
        {
            **config["humas_sd"],
            "input_dir": HUMAS_CENS_SPLIT_DIR,
            "output_dir": HUMAS_SD_OUTDIR,
            "logs_dir": HUMAS_SD_LOGDIR,
            "benchmarks_dir": HUMAS_SD_BMKDIR,
        }


use rule * from HumAS_SD as cens_*


# https://stackoverflow.com/a/63040288
def humas_sd_sm_outputs(wc):
    _ = checkpoints.split_cens_for_humas_sd.get(**wc).output
    wcs = glob_wildcards(
        join(HUMAS_CENS_SPLIT_DIR, wc.sm + "_{chrom}_{ctg}.fa"),
    )
    fnames = [f"{wc.sm}_{chrom}_{ctg}" for chrom, ctg in zip(wcs.chrom, wcs.ctg)]
    chrs = wcs.chrom

    return {
        "stv": expand(rules.cens_generate_stv.output, zip, fname=fnames, chr=chrs),
    }


checkpoint run_humas_sd:
    input:
        rules.cens_generate_monomers.output,
        unpack(humas_sd_sm_outputs),
    output:
        touch(join(HUMAS_SD_OUTDIR, "humas_sd_{sm}.done")),


# https://stackoverflow.com/a/63040288
def humas_sd_chr_outputs(wc):
    _ = [checkpoints.run_humas_sd.get(sm=sm).output for sm in SAMPLE_NAMES]
    wcs = glob_wildcards(
        join(HUMAS_CENS_SPLIT_DIR, "{sm}_{chr}_{ctg}.fa"),
    )
    fnames = [
        f"{sm}_{chrom}_{ctg}"
        for sm, chrom, ctg in zip(wcs.sm, wcs.chr, wcs.ctg)
        if chrom == wc.chr or chrom == f"rc-{wc.chr}"
    ]
    return {
        "stv": [
            *config["plot_hor_stv"].get("ref_stv", []),
            *expand(rules.cens_generate_stv.output, fname=fnames, chr=wc.chr),
        ],
    }


# Force including conda so --containerize includes.
# Must be done since Snakemake won't know rule metadata until runtime.
rule _force_humas_sd_env_inclusion:
    output:
        plots=touch("conda_humas_sd.done"),
    conda:
        "Snakemake-HumAS-SD/workflow/envs/env.yaml"
    shell:
        "echo ''"


rule humas_sd_all:
    input:
        rules._force_humas_sd_env_inclusion.output if IS_CONTAINERIZE_CMD else [],
        expand(rules.run_humas_sd.output, sm=SAMPLE_NAMES),


rule humas_sd_split_cens_only:
    input:
        expand(rules.split_cens_for_humas_sd.output, sm=SAMPLE_NAMES),
