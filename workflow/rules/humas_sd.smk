
include: "common.smk"
include: "utils.smk"


use rule extract_and_index_fa as extract_cens_for_humas_sd with:
    input:
        fa=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}-asm-renamed-reort.fa",
        ),
        bed=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            "interm",
            "{sm}_complete_cens.bed",
        ),
    output:
        seq=temp(
            os.path.join(
                config["humas_sd"]["output_dir"],
                "seq",
                "interm",
                "{sm}_cens.fa",
            )
        ),
        idx=temp(
            os.path.join(
                config["humas_sd"]["output_dir"],
                "seq",
                "interm",
                "{sm}_cens.fa.fai",
            )
        ),
    log:
        "logs/humas_sd/extract_cens_for_humas_sd_{sm}.log",
    params:
        # Seqtk outputs 1-based coords which causes name issues.
        bed=lambda wc, input: f"""<(awk -v OFS="\\t" '{{print $1, $2-1, $3}}' {input.bed})""",
        added_cmds="",


checkpoint split_cens_for_humas_sd:
    input:
        fa=rules.extract_cens_for_humas_sd.output.seq,
    output:
        touch(
            os.path.join(
                config["humas_sd"]["output_dir"],
                "split_cens_for_humas_sd_{sm}.done",
            )
        ),
    log:
        "logs/humas_sd/split_{sm}_cens_for_humas_sd.log",
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
            "logs_dir": "logs/humas_sd",
            "benchmarks_dir": "benchmarks/humas_sd",
        }


use rule * from HumAS_SD as cens_*


# https://stackoverflow.com/a/63040288
def humas_sd_outputs(wc):
    _ = checkpoints.split_cens_for_humas_sd.get(**wc).output
    wcs = glob_wildcards(
        os.path.join(HUMAS_CENS_SPLIT_DIR, wc.sm + "_{chrom}_{ctg}.fa"),
    )
    fnames = [f"{wc.sm}_{chrom}_{ctg}" for chrom, ctg in zip(wcs.chrom, wcs.ctg)]
    chrs = wcs.chrom

    return {
        "hor_bed": expand(
            rules.cens_convert_to_bed9.output, zip, fname=fnames, chr=chrs
        ),
        "stv_row_bed": expand(
            rules.cens_generate_stv.output, zip, fname=fnames, chr=chrs
        ),
    }


checkpoint run_humas_sd:
    input:
        rules.cens_generate_monomers.output,
        unpack(humas_sd_outputs),
    output:
        touch(os.path.join(config["humas_sd"]["output_dir"], "humas_sd_{sm}.done")),


# https://stackoverflow.com/a/63040288
def humas_sd_stv_outputs(wc):
    _ = [checkpoints.run_humas_sd.get(sm=sm).output for sm in SAMPLE_NAMES]
    wcs = glob_wildcards(
        os.path.join(HUMAS_CENS_SPLIT_DIR, "{sm}_" + wc.chr + "_{ctg}.fa"),
    )
    fnames = [f"{sm}_{wc.chr}_{ctg}" for sm, ctg in zip(wcs.sm, wcs.ctg)]
    return {
        "stv": [
            config["plot_hor_stv"]["chm1_stv"],
            config["plot_hor_stv"]["chm13_stv"],
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
