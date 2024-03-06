
rule fmt_correct_alr_regions:
    input:
        alr_fa=rules.extract_correct_alr_regions_rm.output.seq,
        legend=rules.fix_incorrect_merged_legend.output.corrected_legend,
    output:
        fmt_alr_fa=temp(
            os.path.join(
                config["humas_hmmer"]["output_dir"],
                "{sm}_fmt_correct_ALR_regions.500kbp.fa",
            )
        ),
    run:
        with open(str(input.legend)) as legend:
            replacements = dict(tuple(l.strip().split()) for l in legend.readlines())

        with (
            open(str(input.alr_fa)) as input_fa,
            open(str(output.fmt_alr_fa), "wt") as output_fa,
        ):
            for line in input_fa.readlines():
                if line.startswith(">"):
                    record_id, _, _ = line[1:].strip().partition(":")
                    replacement = replacements.get(record_id, record_id)
                    output_fa.write(f">{replacement}\n")
                else:
                    output_fa.write(line)


rule create_rc_merged_legend:
    input:
        rules.fix_incorrect_merged_legend.output.corrected_legend,
    output:
        os.path.join(config["humas_hmmer"]["output_dir"], "{sm}_rc_merged_legend.txt"),
    shell:
        "grep 'rc_chr' {input} > {output}"


use rule fmt_correct_alr_regions as fmt_rc_correct_alr_regions with:
    input:
        alr_fa=rules.rc_correct_alr_regions_rm.output.rc_seq,
        legend=rules.create_rc_merged_legend.output,
    output:
        fmt_alr_fa=temp(
            os.path.join(
                config["humas_hmmer"]["output_dir"],
                "{sm}_fmt_rc_correct_ALR_regions.500kbp.fa",
            )
        ),


rule merge_correct_alr_regions:
    input:
        alr_fa=expand(rules.fmt_correct_alr_regions.output, sm=SAMPLE_NAMES),
    output:
        seq=os.path.join(
            config["humas_hmmer"]["output_dir"], "all_correct_ALR_regions.500kbp.fa"
        ),
        idx=os.path.join(
            config["humas_hmmer"]["output_dir"],
            "all_correct_ALR_regions.500kbp.fa.fai",
        ),
    log:
        "logs/merge_correct_alr_regions.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        cat {input.alr_fa} > {output.seq} 2> {log}
        samtools faidx {output.seq} 2> {log}
        """


use rule merge_correct_alr_regions as merge_correct_alr_regions_rc with:
    input:
        alr_fa=expand(rules.fmt_rc_correct_alr_regions.output, sm=SAMPLE_NAMES),
    output:
        seq=os.path.join(
            config["humas_hmmer"]["output_dir"], "all_correct_ALR_regions.500kbp.rc.fa"
        ),
        idx=os.path.join(
            config["humas_hmmer"]["output_dir"],
            "all_correct_ALR_regions.500kbp.rc.fa.fai",
        ),
    log:
        "logs/merge_rc_correct_alr_regions.log",


rule extract_cens_for_humas_hmmer:
    input:
        all_correct_alr_fa=rules.merge_correct_alr_regions.output.seq,
        rc_all_correct_alr_fa=rules.merge_correct_alr_regions_rc.output.seq,
        corrected_cens_list=rules.split_corrected_rm_output.output.corrected_cens_list,
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
        samtools faidx {output.cens} 2> {log}
        """


checkpoint split_cens_for_humas_hmmer:
    input:
        rules.extract_cens_for_humas_hmmer.output.cens,
    output:
        temp(touch("/tmp/split_cens_for_humas_hmmer_{chr}.done")),
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
        github("koisland/Smk-HumAS-HMMER", path="workflow/Snakefile", branch="main")
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
        temp(touch("/tmp/humas_hmmer_{chr}.done")),
