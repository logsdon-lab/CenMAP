ANNOTATE_REPEATS = config["simplify_repeatmasker_annotations"][
    "repeat_name_pattern_replacement"
]


rule create_annotated_satellites:
    input:
        rules.fix_incorrect_mapped_cens.output.corrected_rm_out,
    output:
        os.path.join(
            config["repeatmasker_simplified"]["output_dir"],
            "all_cens_{repeat}.satellites.bed",
        ),
    params:
        pattern=lambda wc: ANNOTATE_REPEATS[str(wc.repeat)]["pattern"],
        color=lambda wc: ANNOTATE_REPEATS[str(wc.repeat)]["color"],
    shell:
        """
        grep "{params.pattern}" cens_to_annotate.fa.out | \
        sed -e 's/:/\\t/g' -e 's/-/\\t/g' | \
        awk -v OFS="\\t" '{{print $5"-"$6, $7+$9, $7+$10, "{wildcards.repeat}", "0", ".", $7+$9, $7+$10, "{params.color}"}}' > {output}
        """


rule create_ct_track:
    input:
        rules.fix_incorrect_mapped_cens.output.corrected_rm_out,
    output:
        os.path.join(config["repeatmasker_simplified"]["output_dir"], "all_cens.ct.bed"),
    params:
        color=config["simplify_repeatmasker_annotations"]["ct_track_color"],
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
        satellites=rules.create_annotated_satellites.output,
        ct_track=rules.create_ct_track.output,
    output:
        os.path.join(
            config["repeatmasker_simplified"]["output_dir"], "all_cens.annotation.bed"
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
        os.path.join(
            config["repeatmasker_simplified"]["output_dir"], "all_cens.annotation.plot"
        ),
    log:
        "logs/plot_satellite_annotations.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        Rscript {input.script} {input.all_annotations} {output} 2> {log}
        """
