
include: "common.smk"


rule plot_complete_cens:
    input:
        script="workflow/scripts/repeatStructure.R",
        rm_sat_out=os.path.join(
            config["repeatmasker_sat_annot"]["output_dir"],
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
            "all_cens_{chr}_{mer_order}.png",
        ),
        cen_plot_dir=directory(
            os.path.join(
                config["plot_hor_stv"]["output_dir"],
                "plots",
                "{chr}_{mer_order}",
            )
        ),
    log:
        "logs/plot_{chr}_{mer_order}_complete_cens.log",
    conda:
        "../env/r.yaml"
    params:
        hor_filter=0,
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
            --mer_order {wildcards.mer_order} 2> {log}
        fi
        """


rule plot_cens_only:
    input:
        expand(
            rules.plot_complete_cens.output,
            chr=CHROMOSOMES,
            mer_order=config["plot_hor_stv"]["mer_order"],
        ),
