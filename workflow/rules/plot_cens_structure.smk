
include: "common.smk"


rule filter_annotations_cens_structure:
    input:
        all_cdr_bed=lambda wc: (
            os.path.join(
                config["cdr_finder"]["output_dir"],
                "bed",
                "all_cdrs.bed",
            )
            if config.get("cdr_finder")
            else []
        ),
    output:
        cdr_bed=temp(
            os.path.join(
                config["plot_hor_stv"]["output_dir"],
                "bed",
                "{chr}_cdrs.bed",
            )
        ),
    params:
        cdr_output=bool(config.get("cdr_finder", False)),
    log:
        "logs/plot_cens_structure/filter_annotations_{chr}_cens_structure.log",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        if [ "{params.cdr_output}" == "False" ]; then
            touch {output.cdr_bed}
        else
            ( grep '{wildcards.chr}_' {input.all_cdr_bed} || true ) > {output.cdr_bed}
        fi
        """


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
        hor_stv_ort=os.path.join(
            config["plot_hor_stv"]["output_dir"], "bed", "{chr}_AS-HOR_stv_row.ort.bed"
        ),
        cdrs=rules.filter_annotations_cens_structure.output.cdr_bed,
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
            --input_stv_ort {input.hor_stv_ort} \
            --input_cdr {input.cdrs} \
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
    default_target: True
