include: "common.smk"


rule aggregate_all_live_hor:
    input:
        os.path.join(
            config["plot_hor_stv"]["output_dir"], "bed", "results_{chr}_stv"
        ),
    output:
        os.path.join(
            config["cluster_cens"]["output_dir"],
            "bed",
            "{chr}_AS-HOR_liveHORs.all.bed",
        ),
    params:
        liveHORs_pattern=lambda wc, input: os.path.join(
            str(input), "AS-HOR*liveHORs.bed"
        ),
    shell:
        """
        ( cat {params.liveHORs_pattern} || true ) > {output}
        """


rule cluster_cens:
    input:
        script="workflow/scripts/cluster_cens.py",
        live_hor_bed=rules.aggregate_all_live_hor.output,
        cen_img_dir=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "plots",
            f"{{chr}}_{config['cluster_cens']['mer_order']}",
        ),
    output:
        plot=os.path.join(
            config["cluster_cens"]["output_dir"], "{chr}_cen_clusters.png"
        ),
        clusters=os.path.join(
            config["cluster_cens"]["output_dir"], "{chr}_cen_clusters.json"
        ),
    params:
        linkage_method="average",
    log:
        "logs/cluster_cens_{chr}.log",
    conda:
        "../env/py.yaml"
    shell:
        """
        python {input.script} \
        -i {input.live_hor_bed} \
        -p {output.plot} \
        -d {input.cen_img_dir} \
        --output_clusters {output.clusters} \
        --linkage_method {params.linkage_method} 2> {log}
        """


rule cluster_cens_only:
    input:
        expand(rules.cluster_cens.output, chr=CHROMOSOMES),
