
include: "common.smk"


rule extract_cens_for_humas_sd:
    input:
        fa=os.path.join(
            config["repeatmasker"]["output_dir"],
            "seq",
            "all_complete_correct_cens.fa",
        ),
        idx=os.path.join(
            config["repeatmasker"]["output_dir"],
            "seq",
            "all_complete_correct_cens.fa.fai",
        ),
    output:
        cens=os.path.join(config["humas_sd"]["output_dir"], "{chr}_cens.fa"),
        idx=os.path.join(config["humas_sd"]["output_dir"], "{chr}_cens.fa.fai"),
    log:
        "logs/humas_sd/extract_{chr}_cens_for_humas_sd.log",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        seqtk subseq {input.fa} <(grep "{wildcards.chr}[:_]" {input.idx} | cut -f 1) > {output.cens} 2> {log}
        if [ -s "{output.cens}" ]; then
            samtools faidx {output.cens} 2> {log}
        else
            touch {output.idx}
        fi
        """


checkpoint split_cens_for_humas_sd:
    input:
        rules.extract_cens_for_humas_sd.output.cens,
    output:
        touch(
            os.path.join(
                config["humas_sd"]["output_dir"],
                "split_cens_for_humas_sd_{chr}.done",
            )
        ),
    log:
        "logs/humas_sd/split_{chr}_cens_for_humas_sd.log",
    params:
        split_dir=config["humas_sd"]["input_dir"],
    conda:
        "../envs/tools.yaml"
    shell:
        # https://gist.github.com/astatham/621901
        """
        mkdir -p {params.split_dir}
        cat {input} | awk '{{
            if (substr($0, 1, 1)==">") {{
                filename=("{params.split_dir}/" substr($0,2) ".fa")
            }}
            print $0 > filename
        }}' 2> {log}
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
            rules.cens_convert_to_bed9.output, zip, fname=FNAMES, chr=CHRS
        ),
        "stv_row_bed": expand(
            rules.cens_generate_stv.output, zip, fname=FNAMES, chr=CHRS
        ),
    }


checkpoint run_humas_sd:
    input:
        unpack(humas_sd_outputs),
    output:
        touch(os.path.join(config["humas_sd"]["output_dir"], "humas_sd_{chr}.done")),


# https://stackoverflow.com/a/63040288
def humas_sd_stv_outputs(wc):
    _ = checkpoints.run_humas_sd.get(**wc).output
    fnames, chrs = extract_fnames_and_chr(
        rules.cens_filter_hmm_res_overlaps_as_hor.output[0], filter_chr=str(wc.chr)
    )
    return {
        "stv": [
            config["plot_hor_stv"]["chm1_stv"],
            config["plot_hor_stv"]["chm13_stv"],
            *expand(rules.cens_generate_stv.output, zip, fname=FNAMES, chr=CHRS),
        ],
    }


checkpoint create_humas_sd_stv:
    input:
        unpack(humas_sd_stv_outputs),
    output:
        touch(
            os.path.join(
                config["humas_sd"]["output_dir"],
                "create_humas_sd_stv_{chr}.done",
            )
        ),


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
        expand(rules.run_humas_sd.output, chr=CHROMOSOMES),
        expand(rules.create_humas_sd_stv.output, chr=CHROMOSOMES),


rule humas_sd_split_cens_only:
    input:
        expand(rules.split_cens_for_humas_sd.output, chr=CHROMOSOMES),
