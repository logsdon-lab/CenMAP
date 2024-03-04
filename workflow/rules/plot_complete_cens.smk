
rule plot_complete_cens:
    input:
        rm_sat_out=rules.split_rm_satellite_annotations.output,
        hor_stv_out=rules.aggregate_all_stv_row.output,
        chm1_stv=config["plot_hor_stv"]["chm1_stv"],
        chm13_stv=config["plot_hor_stv"]["chm13_stv"],
    output:
        plot=os.path.join("all_cens_{chr}_{mer_order}.png"),
    log:
        "logs/plot_{chr}_{mer_order}_complete_cens.log",
    conda:
        "../env/r.yaml"
    params:
        hor_filter=0,
    shell:
        """
        Rscript \
        --input_rm_sat {input.rm_sat_out} \
        --input_sf {input.hor_stv_out} \
        --input_sf_chm13 {input.chm13_stv} \
        --input_sf_chm1 {input.chm1_stv} \
        --chr {wildcards.chr} \
        --output {output} \
        --hor_filter {params.hor_filter} \
        --mer_order {wildcards.mer_order} 2> {log}
        """
