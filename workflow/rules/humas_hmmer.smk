
rule fmt_correct_alr_regions:
    input:
        alr_fa=rules.extract_correct_alr_regions_rm.output.seq,
        legend=rules.merge_legends_for_rm.output.merged_legends,
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
                    record_id, _, coords = line[1:].strip().partition(":")
                    replacement = replacements.get(record_id, record_id)
                    if coords:
                        output_fa.write(f">{replacement}:{coords}\n")
                    else:
                        output_fa.write(f">{replacement}\n")
                else:
                    output_fa.write(line)


rule create_rc_merged_legend:
    input:
        rules.merge_legends_for_rm.output.merged_legends,
    output:
        os.path.join(config["humas_hmmer"]["output_dir"], "{sm}_rc_merged_legend.txt"),
    shell:
        "sed 's/chr/rc_chr/g' {input} > {output}"


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
        rules.extract_cens_for_humas_hmmer.output.cens,
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
    benchmark:
        "benchmarks/run_humas_hmmer_for_anvil_{chr}.tsv"
    log:
        "logs/run_humas_hmmer_for_anvil_{chr}.log",
    shell:
        """
        mkdir -p {output}
        ./{input.script} {input.input_dir} {output} {input.model} {threads} 2> {log}
        """


rule get_stv_row_from_humas_hmmer_out:
    input:
        humas_hmmer_dir=rules.run_humas_hmmer_for_anvil.output,
        script_filter_live_hor="workflow/scripts/stv/scripts/live_HORs_filter.py",
        script_mon_to_stv="workflow/scripts/stv/scripts/mon2stv.py",
    output:
        directory(
            os.path.join(config["humas_hmmer"]["output_dir"], "results_{chr}_stv")
        ),
    params:
        renamed_bed=lambda wc, output: os.path.join(str(output), "${fname}_renamed.bed"),
        live_hor_bed=lambda wc, output: os.path.join(
            str(output), "${fname}_liveHORs.bed"
        ),
        stv_row_bed=lambda wc, output: os.path.join(str(output), "${fname}_stv_row.bed"),
        as_hor_stv_row_bed=lambda wc, output: os.path.join(
            str(output), "AS-HOR_${fname}_stv_row.bed"
        ),
    log:
        "logs/get_stv_row_from_{chr}_humas_hmmer_out.log",
    conda:
        "../env/py.yaml"
    shell:
        """
        for fname_full in $(find {input.humas_hmmer_dir} -name 'AS-HOR-vs-*.bed'); do
            fname_base=$(basename $fname_full)
            fname="${{fname_base%.bed}}"
            awk -v OFS="\\t" '{{print "{wildcards.chr}", $2, $3, $4, $5, $6, $7, $8, $9, $1}}' $fname_full > "{params.renamed_bed}" 2> {log}
            python3 {input.script_filter_live_hor} "{params.renamed_bed}" > "{params.live_hor_bed}" 2>> {log}
            python3 {input.script_mon_to_stv} "{params.live_hor_bed}" > "{params.stv_row_bed}" 2>> {log}
            awk -v OFS="\\t" 'FNR==NR{{a[NR]=$10;next}}{{$1=a[FNR]}}1' "{params.renamed_bed}" "{params.stv_row_bed}" > "{params.as_hor_stv_row_bed}" 2>> {log}
        done
        """


rule aggregate_all_stv_row:
    input:
        rules.get_stv_row_from_humas_hmmer_out.output,
    output:
        os.path.join(
            config["humas_hmmer"]["output_dir"], "{chr}_AS-HOR_stv_row.all.bed"
        ),
    params:
        stv_row_pattern=lambda wc, input: os.path.join(
            str(input), "AS-HOR_*_stv_row.bed"
        ),
    shell:
        """
        cat {params.stv_row_pattern} > {output}
        """
