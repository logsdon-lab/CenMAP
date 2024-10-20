
include: "common.smk"


rule plot_cens_structure:
    input:
        script="workflow/scripts/plot_cens_all_stv.R",
        rm_sat_out=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "repeats",
            "all_cens_{chr}.annotation.fa.out",
        ),
        hor_stv_out=os.path.join(
            config["plot_hor_stv"]["output_dir"], "bed", "{chr}_AS-HOR_stv_row.all.bed"
        ),
        chm1_stv=config["plot_hor_stv"]["chm1_stv"],
        chm13_stv=config["plot_hor_stv"]["chm13_stv"],
    output:
        plot=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "plots",
            "hor",
            "all_cens_{chr}.png",
        ),
        cen_plot_dir=directory(
            os.path.join(
                config["plot_hor_stv"]["output_dir"],
                "plots",
                "hor",
                "{chr}",
            )
        ),
    log:
        "logs/plot_cens_structure/plot_{chr}_cens_structure.log",
    conda:
        "../envs/r.yaml"
    params:
        hor_filter=0,
        mer_order=lambda wc: MONOMER_ORDER[wc.chr],
    shell:
        """
        if ! [ -s  {input.hor_stv_out} ] || ! [ -s {input.rm_sat_out} ]; then
            touch {output.plot}
            mkdir -p {output.cen_plot_dir}
        else
            Rscript {input.script} \
            --input_rm_sat {input.rm_sat_out} \
            --input_stv {input.hor_stv_out} \
            --input_stv_chm13 {input.chm13_stv} \
            --input_stv_chm1 {input.chm1_stv} \
            --chr {wildcards.chr} \
            --output {output.plot} \
            --output_dir {output.cen_plot_dir} \
            --hor_filter {params.hor_filter} \
            --mer_order {params.mer_order} 2> {log}
        fi
        """


rule plot_cens_structure_only:
    input:
        expand(
            rules.plot_cens_structure.output,
            chr=CHROMOSOMES,
        ),
