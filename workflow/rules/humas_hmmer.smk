
include: "common.smk"


rule extract_cens_for_humas_hmmer:
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
        cens=os.path.join(config["humas_hmmer"]["output_dir"], "{chr}_cens.fa"),
        idx=os.path.join(config["humas_hmmer"]["output_dir"], "{chr}_cens.fa.fai"),
    log:
        "logs/humas_hmmer/extract_{chr}_cens_for_humas_hmmer.log",
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


checkpoint split_cens_for_humas_hmmer:
    input:
        rules.extract_cens_for_humas_hmmer.output.cens,
    output:
        touch(
            os.path.join(
                config["humas_hmmer"]["output_dir"],
                "split_cens_for_humas_hmmer_{chr}.done",
            )
        ),
    log:
        "logs/humas_hmmer/split_{chr}_cens_for_humas_hmmer.log",
    params:
        split_dir=config["humas_hmmer"]["input_dir"],
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


MODULE_REPO = "logsdon-lab/Snakemake-HumAS-HMMER"
MODULE_BRANCH = "main"
MODULE_PATH = "workflow/Snakefile"


module HumAS_HMMER:
    snakefile:
        github(
            MODULE_REPO,
            path=MODULE_PATH,
            branch=MODULE_BRANCH,
        )
    config:
        config["humas_hmmer"]


use rule * from HumAS_HMMER as cens_*


MODULE_ENV_FILE = rules.cens_humas_hmmer_analysis.__dict__["rule"]._conda_env


# https://stackoverflow.com/a/63040288
def humas_hmmer_outputs(wc):
    _ = checkpoints.split_cens_for_humas_hmmer.get(**wc).output
    fnames, chrs = extract_fnames_and_chr(
        os.path.join(config["humas_hmmer"]["input_dir"], "{fname}.fa"),
        filter_chr=str(wc.chr),
    )
    return {
        "overlaps": expand(
            rules.cens_filter_hmm_res_overlaps_as_hor.output, fname=fnames
        ),
        "final": expand(rules.cens_filter_final_hmm_res_as_hor.output, fname=fnames),
    }


checkpoint run_humas_hmmer_for_anvil:
    input:
        unpack(humas_hmmer_outputs),
    output:
        touch(
            os.path.join(config["humas_hmmer"]["output_dir"], "humas_hmmer_{chr}.done")
        ),


rule get_live_hor:
    input:
        script="workflow/scripts/stv_fix/scripts/live_HORs_filter.py",
        # fname is wildcard
        humas_hmmer_out=rules.cens_filter_hmm_res_overlaps_as_hor.output,
    output:
        renamed_bed=os.path.join(
            config["humas_hmmer"]["output_dir"],
            "bed",
            "{chr}",
            "{fname}_renamed.bed",
        ),
        live_hor_bed=os.path.join(
            config["humas_hmmer"]["output_dir"],
            "bed",
            "{chr}",
            "{fname}_liveHORs.bed",
        ),
    params:
        # Fix inversion for chromosome 1 and 19.
        inversion_fix_cmd=lambda wc: (
            "| sed 's+S1C1/5/19H1L.6/4+S1C1/5/19H1L.6+g'"
            if wc.chr == "chr1" or wc.chr == "chr19"
            else ""
        ),
    log:
        "logs/humas_hmmer/get_live_hor_{fname}_{chr}.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        awk -v OFS="\\t" '{{
            print "{wildcards.chr}", $2, $3, $4, $5, $6, $7, $8, $9, $1
        }}' {input.humas_hmmer_out} > {output.renamed_bed} 2> {log}
        {{ python3 {input.script} {output.renamed_bed} {params.inversion_fix_cmd} ;}} > {output.live_hor_bed} 2>> {log}
        """


rule filter_as_hor_stv_bed:
    input:
        script="workflow/scripts/stv_fix/scripts/mon2stv.py",
        live_hor_bed=rules.get_live_hor.output.live_hor_bed,
        renamed_bed=rules.get_live_hor.output.renamed_bed,
    output:
        stv_row_bed=os.path.join(
            config["humas_hmmer"]["output_dir"],
            "bed",
            "{chr}",
            "{fname}_stv_row.bed",
        ),
        as_hor_stv_row_bed=os.path.join(
            config["humas_hmmer"]["output_dir"],
            "bed",
            "{chr}",
            "AS-HOR_{fname}_stv_row.bed",
        ),
    log:
        "logs/humas_hmmer/filter_as_hor_stv_bed_{fname}_{chr}.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        python3 {input.script} {input.live_hor_bed} > {output.stv_row_bed} 2>> {log}
        awk -v OFS="\\t" 'FNR==NR{{a[NR]=$10;next}}{{$1=a[FNR]}}1' {input.renamed_bed} {output.stv_row_bed} > {output.as_hor_stv_row_bed} 2>> {log}
        """


# https://stackoverflow.com/a/63040288
def humas_hmmer_stv_outputs(wc):
    _ = checkpoints.run_humas_hmmer_for_anvil.get(**wc).output
    fnames, chrs = extract_fnames_and_chr(
        rules.cens_filter_hmm_res_overlaps_as_hor.output[0], filter_chr=str(wc.chr)
    )
    return {
        "stv": [
            config["plot_hor_stv"]["chm1_stv"],
            config["plot_hor_stv"]["chm13_stv"],
            *expand(
                rules.filter_as_hor_stv_bed.output.as_hor_stv_row_bed,
                zip,
                fname=fnames,
                chr=chrs,
            ),
        ],
    }


checkpoint create_humas_hmmer_stv:
    input:
        unpack(humas_hmmer_stv_outputs),
    output:
        touch(
            os.path.join(
                config["humas_hmmer"]["output_dir"],
                "create_humas_hmmer_stv_{chr}.done",
            )
        ),


# Force including conda so --containerize includes.
# Must be done since Snakemake won't know rule metadata until runtime.
rule _force_humas_hmmer_env_inclusion:
    output:
        plots=touch("conda_humas_hmmer.done"),
    conda:
        f"https://raw.githubusercontent.com/{MODULE_REPO}/{MODULE_BRANCH}/workflow/{MODULE_ENV_FILE}"
    shell:
        "echo ''"


rule humas_hmmer_only:
    input:
        rules._force_humas_hmmer_env_inclusion.output if IS_CONTAINERIZE_CMD else [],
        expand(rules.run_humas_hmmer_for_anvil.output, chr=CHROMOSOMES),
        expand(rules.create_humas_hmmer_stv.output, chr=CHROMOSOMES),


rule humas_hmmer_split_cens_only:
    input:
        expand(rules.split_cens_for_humas_hmmer.output, chr=CHROMOSOMES),
