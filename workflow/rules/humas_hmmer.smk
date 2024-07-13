
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
            "all_correct_cens.fa.fai",
        ),
    output:
        cens=os.path.join(config["humas_hmmer"]["output_dir"], "{chr}_cens.fa"),
        idx=os.path.join(config["humas_hmmer"]["output_dir"], "{chr}_cens.fa.fai"),
    log:
        "logs/humas_hmmer/extract_{chr}_cens_for_humas_hmmer.log",
    conda:
        "../env/tools.yaml"
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
        "../env/tools.yaml"
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


module HumAS_HMMER:
    snakefile:
        github(
            "logsdon-lab/Snakemake-HumAS-HMMER",
            path="workflow/Snakefile",
            branch="main",
        )
    config:
        config["humas_hmmer"]


use rule * from HumAS_HMMER as cens_*


# https://stackoverflow.com/a/63040288
def humas_hmmer_outputs(wc):
    _ = checkpoints.split_cens_for_humas_hmmer.get(**wc).output
    fnames = glob_wildcards(
        os.path.join(config["humas_hmmer"]["input_dir"], "{fname}.fa")
    ).fname
    return {
        "overlaps": expand(
            rules.cens_filter_hmm_res_overlaps_as_hor.output, fname=fnames
        ),
        "final": expand(rules.cens_filter_final_hmm_res_as_hor.output, fname=fnames),
    }


rule run_humas_hmmer_for_anvil:
    input:
        unpack(humas_hmmer_outputs),
    output:
        touch(
            os.path.join(config["humas_hmmer"]["output_dir"], "humas_hmmer_{chr}.done")
        ),


rule humas_hmmer_only:
    input:
        expand(rules.run_humas_hmmer_for_anvil.output, chr=CHROMOSOMES),
