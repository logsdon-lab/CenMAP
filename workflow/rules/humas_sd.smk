
include: "common.smk"


# TODO: Needs to be after repeatmasker correction.
checkpoint split_cens_for_humas_sd:
    input:
        fa=lambda wc: expand(
            rules.extract_alr_region_sample_by_chr.output.seq,
            sm=SAMPLE_NAMES,
            chr=wc.chr,
        ),
        rename_key=lambda wc: expand(
            rules.map_collapse_cens.output.renamed_cens_key, sm=SAMPLE_NAMES
        ),
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
        awk '{{
            # Read key values in first file.
            if (FNR == NR) {{
                # Add coords to name.
                kv[$1]=$2;
                next;
            }}
            if (substr($0, 1, 1)==">") {{
                ctg_name=substr($0,2)
                split(ctg_name, ctg_name_parts, ":")
                new_ctg_name=kv[ctg_name_parts[1]]":"ctg_name_parts[2]
                filename=("{params.split_dir}/" new_ctg_name ".fa")
            }}
            print $0 > filename
        }}' <(awk -v OFS="\\t" '$5=="{wildcards.chr}"' {input.rename_key}) <(cat {input.fa}) 2> {log}
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
        touch(os.path.join(config["humas_sd"]["output_dir"], "humas_sd_{chr}.done")),


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
