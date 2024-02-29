ANNOTATE_SAT_REPEATS = config["repeatmasker_sat_annot"][
    "repeat_name_pattern_replacement"
]


rule create_annotated_satellites:
    input:
        rules.fix_incorrect_mapped_cens.output.corrected_rm_out,
    output:
        os.path.join(
            config["repeatmasker_sat_annot"]["output_dir"],
            "all_cens_{repeat}.satellites.bed",
        ),
    params:
        pattern=lambda wc: ANNOTATE_SAT_REPEATS[str(wc.repeat)]["pattern"],
        color=lambda wc: ANNOTATE_SAT_REPEATS[str(wc.repeat)]["color"],
    shell:
        """
        grep "{params.pattern}" {input} | \
        sed -e 's/:/\\t/g' -e 's/-/\\t/g' | \
        awk -v OFS="\\t" '{{print $5"-"$6, $7+$9, $7+$10, "{wildcards.repeat}", "0", ".", $7+$9, $7+$10, "{params.color}"}}' > {output}
        """


rule create_ct_track:
    input:
        rules.fix_incorrect_mapped_cens.output.corrected_rm_out,
    output:
        os.path.join(config["repeatmasker_sat_annot"]["output_dir"], "all_cens.ct.bed"),
    params:
        color=config["repeatmasker_sat_annot"]["ct_track_color"],
    log:
        "logs/create_ct_track.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        sed -e 's/:/\\t/g' -e 's/-/\\t/g' {input} | \
        awk -v OFS="\\t" '{{print $5"-"$6, $7, $8, "ct", "0", ".", $7, $8, "{params.color}"}}' | \
        sort | uniq | grep "chr" > {output}
        """


rule aggregate_rm_satellite_annotations:
    input:
        satellites=expand(
            rules.create_annotated_satellites.output,
            repeat=ANNOTATE_SAT_REPEATS.keys(),
        ),
        ct_track=rules.create_ct_track.output,
    output:
        os.path.join(
            config["repeatmasker_sat_annot"]["output_dir"], "all_cens.annotation.bed"
        ),
    shell:
        """
        cat {input.ct_track} {input.satellites} > {output}
        """


rule plot_satellite_annotations:
    input:
        script="workflow/scripts/repeatStructure_satellites.R",
        all_annotations=rules.aggregate_rm_satellite_annotations.output,
    output:
        chr_annot=temp(
            os.path.join(
                config["repeatmasker_sat_annot"]["output_dir"],
                "all_cens_{chr}.annotation.fa.out",
            )
        ),
        chr_plot=os.path.join(
            config["repeatmasker_sat_annot"]["output_dir"],
            "all_cens_{chr}.annotation.png",
        ),
    log:
        "logs/plot_{chr}_satellite_annotations.log",
    conda:
        "../env/r.yaml"
    shell:
        """
        grep "{wildcards.chr}[_:]" {input.all_annotations} > {output.chr_annot}
        Rscript {input.script} {output.chr_annot} {output.chr_plot} 2> {log}
        """
