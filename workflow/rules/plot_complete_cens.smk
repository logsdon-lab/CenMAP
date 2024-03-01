
rule plot_complete_cen:
    input:
        rm_out=rules.split_corrected_rm_output.output.corrected_rm_out,
        rm_sat_out=rules.split_rm_satellite_annotations.output,
        hor_stv_out=rules.aggregate_all_stv_row.output,
    output:
        plot=os.path.join("all_cens_{chr}.png"),
    log:
        "logs/plot_{chr}_complete_cens.log",
    conda:
        "../env/r.yaml"
    shell:
        """
        """
