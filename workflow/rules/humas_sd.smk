
include: "common.smk"


use rule extract_and_index_fa as extract_cens_for_humas_sd with:
    input:
        fa=os.path.join(
            config["humas_sd"]["output_dir"],
            "{sm}",
            "{sm}_regions.renamed.reort.fa",
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


checkpoint split_cens_for_humas_sd:
    input:
        fa=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}",
            "{sm}_regions.renamed.reort.fa",
        ),
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
        split_dir=config["humas_sd"]["input_dir"],
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
        config["humas_sd"]


use rule * from HumAS_SD as cens_*


# https://stackoverflow.com/a/63040288
def humas_sd_outputs(wc):
    _ = checkpoints.split_cens_for_humas_sd.get(**wc).output
    fnames, chrs = extract_fnames_and_chr(
        os.path.join(config["humas_sd"]["input_dir"], "{fname}.fa"),
        filter_chr=str(wc.chr),
    )
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
        touch(
            os.path.join(config["humas_sd"]["output_dir"], "humas_sd_{sm}_{chr}.done")
        ),


# https://stackoverflow.com/a/63040288
def humas_sd_stv_outputs(wc):
    _ = checkpoints.run_humas_sd.get(**wc).output
    fnames, chrs = extract_fnames_and_chr(
        os.path.join(config["humas_sd"]["input_dir"], "{fname}.fa"),
        filter_chr=str(wc.chr),
    )
    return {
        "stv": [
            config["plot_hor_stv"]["chm1_stv"],
            config["plot_hor_stv"]["chm13_stv"],
            *expand(rules.cens_generate_stv.output, zip, fname=fnames, chr=chrs),
        ],
    }


checkpoint create_humas_sd_stv:
    input:
        unpack(humas_sd_stv_outputs),
    output:
        touch(
            os.path.join(
                config["humas_sd"]["output_dir"],
                "create_humas_sd_stv_{sm}_{chr}.done",
            )
        ),


def humas_sd_stv_sm_outputs(wc):
    _ = [
        checkpoints.run_humas_sd.get(**{"sm": wc.sm, "chr": chrom}).output
        for chrom in CHROMOSOMES
    ]
    wcs = glob_wildcards(
        os.path.join(config["humas_sd"]["input_dir"], "{sm}_{chr}_{ctg_name}.fa")
    )
    outputs = []
    for sm, chrom, ctg_name in zip(wcs.sm, wcs.chr, wcs.ctg_name):
        if sm != wc.sm:
            continue

        outputs.extend(
            expand(
                rules.cens_generate_stv.output,
                zip,
                fname=f"{sm}_{chrom}_{ctg_name}",
                chr=[chrom],
            )
        )
    return outputs


# Force including conda so --containerize includes.
# Must be done since Snakemake won't know rule metadata until runtime.
rule _force_humas_sd_env_inclusion:
    output:
        plots=touch("conda_humas_sd.done"),
    conda:
        "workflow/rules/Snakemake-HumAS-SD/workflow/envs/env.yaml"
    shell:
        "echo ''"


rule humas_sd_all:
    input:
        rules._force_humas_sd_env_inclusion.output if IS_CONTAINERIZE_CMD else [],
        expand(rules.run_humas_sd.output, sm=SAMPLE_NAMES, chr=CHROMOSOMES),
        expand(rules.create_humas_sd_stv.output, sm=SAMPLE_NAMES, chr=CHROMOSOMES),


rule humas_sd_split_cens_only:
    input:
        expand(rules.split_cens_for_humas_sd.output, sm=SAMPLE_NAMES, chr=CHROMOSOMES),
