
include: "common.smk"


rule fmt_correct_alr_regions:
    input:
        alr_fa=os.path.join(
            config["repeatmasker"]["output_dir"],
            "seq",
            "{sm}_correct_ALR_regions.500kbp.fa",
        ),
        legend=os.path.join(
            config["repeatmasker"]["output_dir"],
            "status",
            "{sm}_corrected_merged_legend.txt",
        ),
    output:
        fmt_alr_fa=temp(
            os.path.join(
                config["humas_hmmer"]["output_dir"],
                "{sm}_fmt_correct_ALR_regions.{comp}.500kbp.fa",
            )
        ),
    run:
        with open(str(input.legend)) as legend:
            replacements = dict(tuple(l.strip().split()) for l in legend.readlines())

        with (
            open(str(input.alr_fa)) as input_fa,
            open(str(output.fmt_alr_fa), "wt") as output_fa,
        ):
            skip_record = False
            for line in input_fa.readlines():
                if line.startswith(">"):
                    # HG00171_chr22_h1tg000027l#1-26260313
                    sm, chrm, record_id = line[1:].strip().rsplit("_", maxsplit=2)
                    record_id, *delim_coords = record_id.partition(":")
                    replacement = replacements.get(record_id, "")
                    is_rc = "rc" in replacement
                    if (
                        (wildcards.comp == "rc" and not is_rc)
                        or (wildcards.comp == "c" and is_rc)
                        or replacement == ""
                    ):
                        skip_record = True
                        continue
                    output_fa.write(f">{replacement}{''.join(delim_coords)}\n")
                else:
                    if skip_record:
                        skip_record = False
                        continue
                    output_fa.write(line)


rule merge_correct_alr_regions:
    input:
        alr_fa=lambda wc: expand(
            rules.fmt_correct_alr_regions.output, sm=SAMPLE_NAMES, comp=[wc.comp]
        ),
    output:
        seq=os.path.join(
            config["humas_hmmer"]["output_dir"],
            "all_correct_ALR_regions.{comp}.500kbp.fa",
        ),
        idx=os.path.join(
            config["humas_hmmer"]["output_dir"],
            "all_correct_ALR_regions.{comp}.500kbp.fa.fai",
        ),
    log:
        "logs/merge_correct_alr_regions_{comp}.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        cat {input.alr_fa} > {output.seq} 2> {log}
        if [-s {output.seq}]; then
            samtools faidx {output.seq} 2> {log}
        else
            touch {output.idx}
        fi
        """


rule extract_cens_for_humas_hmmer:
    input:
        all_correct_alr_fa=expand(
            rules.merge_correct_alr_regions.output.seq, comp=["c"]
        ),
        rc_all_correct_alr_fa=expand(
            rules.merge_correct_alr_regions.output.seq, comp=["rc"]
        ),
        corrected_cens_list=os.path.join(
            config["repeatmasker"]["output_dir"], "status", "corrected_{chr}_cens.list"
        ),
    output:
        cens=os.path.join(config["humas_hmmer"]["output_dir"], "{chr}_cens.fa"),
        idx=os.path.join(config["humas_hmmer"]["output_dir"], "{chr}_cens.fa.fai"),
    log:
        "logs/extract_{chr}_cens_for_humas_hmmer.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        seqtk subseq {input.all_correct_alr_fa} {input.corrected_cens_list} > {output.cens} 2> {log}
        seqtk subseq {input.rc_all_correct_alr_fa} {input.corrected_cens_list} >> {output.cens} 2> {log}
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
        "logs/split_{chr}_cens_for_humas_hmmer.log",
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
        # TODO: Sort maybe?
        temp(
            touch(
                os.path.join(
                    config["humas_hmmer"]["input_dir"], "humas_hmmer_{chr}.done"
                )
            )
        ),


rule humas_hmmer_only:
    input:
        expand(rules.fmt_correct_alr_regions.output, sm=SAMPLE_NAMES, comp=["c", "rc"]),
        expand(rules.merge_correct_alr_regions.output, comp=["c", "rc"]),
        expand(rules.extract_cens_for_humas_hmmer.output, chr=CHROMOSOMES),
        expand(rules.split_cens_for_humas_hmmer.output, chr=CHROMOSOMES),
        expand(rules.run_humas_hmmer_for_anvil.output, chr=CHROMOSOMES),
