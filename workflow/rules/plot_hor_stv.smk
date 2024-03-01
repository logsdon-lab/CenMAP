
rule get_stv_row_from_humas_hmmer_out:
    input:
        humas_hmmer_done=rules.run_humas_hmmer_for_anvil.output,
        script_filter_live_hor="workflow/scripts/stv_fix/scripts/live_HORs_filter.py",
        script_mon_to_stv="workflow/scripts/stv_fix/scripts/mon2stv.py",
    output:
        directory(
            os.path.join(config["plot_hor_stv"]["output_dir"], "results_{chr}_stv")
        ),
    params:
        humas_hmmer_dir=lambda wc: config["humas_hmmer"]["output_dir"],
        renamed_bed=lambda wc, output: os.path.join(str(output), "${fname}_renamed.bed"),
        live_hor_bed=lambda wc, output: os.path.join(
            str(output), "${fname}_liveHORs.bed"
        ),
        stv_row_bed=lambda wc, output: os.path.join(str(output), "${fname}_stv_row.bed"),
        as_hor_stv_row_bed=lambda wc, output: os.path.join(
            str(output), "AS-HOR_${fname}_stv_row.bed"
        ),
        # Fix inversion for chromosome 1 and 19.
        inversion_fix_cmd=lambda wc: (
            "| sed 's+S1C1/5/19H1L.6/4+S1C1/5/19H1L.6+g'"
            if wc.chr == "chr1" or wc.chr == "chr19"
            else ""
        ),
    log:
        "logs/get_stv_row_from_{chr}_humas_hmmer_out.log",
    conda:
        "../env/py.yaml"
    shell:
        """
        mkdir -p {output}
        for fname_full in $(find {params.humas_hmmer_dir} -name 'AS-HOR-vs-*{wildcards.chr}_*.bed'); do
            fname_base=$(basename $fname_full)
            fname="${{fname_base%.bed}}"
            awk -v OFS="\\t" '{{print "{wildcards.chr}", $2, $3, $4, $5, $6, $7, $8, $9, $1}}' $fname_full > "{params.renamed_bed}" 2> {log}
            {{ python3 {input.script_filter_live_hor} "{params.renamed_bed}" {params.inversion_fix_cmd} ;}} > "{params.live_hor_bed}" 2>> {log}
            python3 {input.script_mon_to_stv} "{params.live_hor_bed}" > "{params.stv_row_bed}" 2>> {log}
            awk -v OFS="\\t" 'FNR==NR{{a[NR]=$10;next}}{{$1=a[FNR]}}1' "{params.renamed_bed}" "{params.stv_row_bed}" > "{params.as_hor_stv_row_bed}" 2>> {log}
        done
        """


rule aggregate_all_stv_row:
    input:
        rules.get_stv_row_from_humas_hmmer_out.output,
    output:
        os.path.join(
            config["plot_hor_stv"]["output_dir"], "{chr}_AS-HOR_stv_row.all.bed"
        ),
    params:
        stv_row_pattern=lambda wc, input: os.path.join(
            str(input), "AS-HOR_*_stv_row.bed"
        ),
    shell:
        """
        cat {params.stv_row_pattern} > {output}
        """


rule plot_stv_with_order:
    input:
        script="workflow/scripts/StVHORorganization_simplified_2datasets_basedOnStartStop_newcolors.R",
        chm1_stv=config["plot_hor_stv"]["chm1_stv"],
        chm13_stv=config["plot_hor_stv"]["chm13_stv"],
        all_stv=rules.aggregate_all_stv_row.output,
    output:
        hor_array_plot=os.path.join(
            config["plot_hor_stv"]["output_dir"], "{chr}_hgsvc3_{mer_order}ontop.png"
        ),
        hor_array_length=os.path.join(
            config["plot_hor_stv"]["output_dir"], "{chr}_length_{mer_order}.tsv"
        ),
    log:
        "logs/plot_{chr}_stv_{mer_order}_on_top.log",
    conda:
        "../env/r.yaml"
    shell:
        """
        Rscript {input.script} \
        --input {input.all_stv} \
        --input_chm13 {input.chm13_stv} \
        --input_chm1 {input.chm1_stv} \
        --output {output.hor_array_plot} \
        --output_tbl {output.hor_array_length} \
        --chr {wildcards.chr} \
        --mer_order {wildcards.mer_order} 2> {log}
        """
