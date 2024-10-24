include: "common.smk"


rule merge_complete_and_correct_rm_out:
    input:
        expand(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "repeats",
                "all",
                "complete_correct_{chr}_cens.fa.out",
            ),
            chr=CHROMOSOMES,
        ),
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "all_samples_and_ref_complete_correct_cens.fa.out",
        ),
    shell:
        """
        cat {input} > {output}
        """


rule create_annotated_satellites:
    input:
        ref_rm_out=config["repeatmasker"]["ref_repeatmasker_output"],
        corrected_rm_out=rules.merge_complete_and_correct_rm_out.output,
    output:
        os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "bed",
            "all_cens_{repeat}.satellites.bed",
        ),
    conda:
        "../envs/tools.yaml"
    params:
        pattern=lambda wc: ANNOTATE_SAT_REPEATS[str(wc.repeat)]["pattern"],
        color=lambda wc: ANNOTATE_SAT_REPEATS[str(wc.repeat)]["color"],
    log:
        "logs/plot_repeatmasker_sat_annot/create_{repeat}_annotated_satellites.log",
    shell:
        """
        {{ cat {input} | \
        grep "{params.pattern}" | \
        awk -v OFS="\\t" '{{
            # Find start in contig name.
            match($5, ":(.+)-", starts);
            # Create adjust start/end
            new_start=$6+starts[1];
            new_end=$7+starts[1];

            print $5, new_start, new_end, "{wildcards.repeat}", "0", ".", new_start, new_end, "{params.color}"
        }}';}} > {output} 2> {log}
        """


rule create_ct_track:
    input:
        ref_rm_out=config["repeatmasker"]["ref_repeatmasker_output"],
        corrected_rm_out=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "all_samples_and_ref_complete_correct_cens.fa.out",
        ),
    output:
        os.path.join(config["plot_hor_stv"]["output_dir"], "bed", "all_cens.ct.bed"),
    params:
        color=ANNOTATE_SAT_REPEATS["ct"]["color"],
    log:
        "logs/plot_repeatmasker_sat_annot/create_ct_track.log",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ cat {input} | \
        awk -v OFS="\\t" '{{
            # Find start in contig name.
            match($5, ":(.+)-", starts);
            match($5, ".*-(.+)$", ends);
            print $5, starts[1], ends[1], "ct", "0", ".", starts[1], ends[1], "{params.color}"
        }}' | \
        sort | uniq | grep -P "chr|cen";}} > {output} 2> {log}
        """


rule aggregate_rm_satellite_annotations:
    input:
        satellites=expand(
            rules.create_annotated_satellites.output,
            repeat=[r for r in ANNOTATE_SAT_REPEATS.keys() if r != "ct"],
        ),
        ct_track=rules.create_ct_track.output,
    output:
        os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "bed",
            "all_cens.annotation.bed",
        ),
    shell:
        """
        cat {input.ct_track} {input.satellites} > {output}
        """


rule split_rm_satellite_annotations:
    input:
        all_annotations=rules.aggregate_rm_satellite_annotations.output,
    output:
        chr_annot=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "repeats",
            "all_cens_{chr}.annotation.fa.out",
        ),
    params:
        chr_pattern=lambda wc: str(wc.chr).replace("chr", r"(chr|cen)") + r"[_:\-\\tv]",
    shell:
        """
        grep -P "{params.chr_pattern}" {input.all_annotations} > {output.chr_annot}
        """


rule plot_satellite_annotations:
    input:
        script="workflow/scripts/plot_cens_onlysatRM.R",
        chr_annot=rules.split_rm_satellite_annotations.output,
    output:
        chr_plot=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "plots",
            "satellite",
            "all_cens_{chr}.annotation.png",
        ),
    log:
        "logs/plot_repeatmasker_sat_annot/plot_{chr}_satellite_annotations.log",
    conda:
        "../envs/r.yaml"
    shell:
        """
        Rscript {input.script} {input.chr_annot} {output.chr_plot} 2> {log}
        """


rule plot_repeatmasker_sat_only:
    input:
        expand(rules.plot_satellite_annotations.output, chr=CHROMOSOMES),
    default_target: True
