
rule extract_cens_for_humas_hmmer:
    input:
        all_correct_alr_fa=rules.format_add_censat_annot_repeatmasker_output.output,
        rc_all_correct_alr_fa=rules.reverse_complete_repeatmasker_output.output,
        corrected_cens_list=rules.create_correct_oriented_cens_list.output.corrected_cens_list,
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


rule split_cens_for_humas_hmmer:
    input:
        rules.extract_cens_for_humas_hmmer.output,
    output:
        directory(os.path.join(config["humas_hmmer"]["output_dir"], "{chr}")),
    log:
        "logs/split_{chr}_cens_for_humas_hmmer.log",
    conda:
        "../env/tools.yaml"
    shell:
        # https://gist.github.com/astatham/621901
        """
        mkdir -p {output}
        cat {input} | awk '{{
            if (substr($0, 1, 1)==">") {{
                filename=("{output}/" substr($0,2) ".fa")
            }}
            print $0 > filename
        }}' 2> {log}
        """


rule run_humas_hmmer_for_anvil:
    input:
        script="workflow/scripts/HumAS-HMMER_for_AnVIL/hmmer-run.sh",
        input_dir=rules.split_cens_for_humas_hmmer.output,
        model=config["humas_hmmer"]["model"],
    output:
        directory(os.path.join(config["humas_hmmer"]["output_dir"], "results_{chr}")),
    conda:
        "../env/tools.yaml"
    threads: config["humas_hmmer"]["threads"]
    log:
        "logs/run_humas_hmmer_for_anvil_{chr}.log",
    shell:
        """
        ./{input.script} {input.input_dir} {output} {input.model}
        """
