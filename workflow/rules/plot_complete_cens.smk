
include: "common.smk"


rule plot_complete_cens:
    input:
        script="workflow/scripts/repeatStructure.R",
        rm_sat_out=os.path.join(
            config["repeatmasker_sat_annot"]["output_dir"],
            "all_cens_{chr}.annotation.fa.out",
        ),
        hor_stv_out=os.path.join(
            config["plot_hor_stv"]["output_dir"], "{chr}_AS-HOR_stv_row.all.bed"
        ),
        chm1_stv=config["plot_hor_stv"]["chm1_stv"],
        chm13_stv=config["plot_hor_stv"]["chm13_stv"],
    output:
        plot=os.path.join(
            config["plot_hor_stv"]["output_dir"], "all_cens_{chr}_{mer_order}.png"
        ),
    log:
        "logs/plot_{chr}_{mer_order}_complete_cens.log",
    conda:
        "../env/r.yaml"
    params:
        hor_filter=0,
    shell:
        """
        Rscript {input.script} \
        --input_rm_sat {input.rm_sat_out} \
        --input_sf {input.hor_stv_out} \
        --input_sf_chm13 {input.chm13_stv} \
        --input_sf_chm1 {input.chm1_stv} \
        --chr {wildcards.chr} \
        --output {output} \
        --hor_filter {params.hor_filter} \
        --mer_order {wildcards.mer_order} 2> {log}
        """


rule plot_cens_only:
    input:
        expand(
            rules.plot_complete_cens.output,
            chr=CHROMOSOMES,
            mer_order=config["plot_hor_stv"]["mer_order"],
        ),
