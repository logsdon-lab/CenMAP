

include: "common.smk"


ANNOTATE_SAT_REPEATS = config["repeatmasker_sat_annot"][
    "repeat_name_pattern_replacement"
]


rule create_annotated_satellites:
    input:
        ref_rm_out=(
            config["repeatmasker"]["ref_repeatmasker_output"]
            if config["repeatmasker"].get("ref_repeatmasker_output")
            else rules.run_repeatmasker_ref.output
        ),
        # corrected_rm_out=rules.fix_incorrect_mapped_cens.output.corrected_rm_out,
        corrected_rm_out=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "all_corrected_cens.fa.out",
        ),
    output:
        os.path.join(
            config["repeatmasker_sat_annot"]["output_dir"],
            "bed",
            "all_cens_{repeat}.satellites.bed",
        ),
    params:
        pattern=lambda wc: ANNOTATE_SAT_REPEATS[str(wc.repeat)]["pattern"],
        color=lambda wc: ANNOTATE_SAT_REPEATS[str(wc.repeat)]["color"],
    log:
        "logs/create_{repeat}_annotated_satellites.log",
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
        ref_rm_out=(
            config["repeatmasker"]["ref_repeatmasker_output"]
            if config["repeatmasker"].get("ref_repeatmasker_output")
            else rules.run_repeatmasker_ref.output
        ),
        corrected_rm_out=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "all_corrected_cens.fa.out",
        ),
    output:
        os.path.join(
            config["repeatmasker_sat_annot"]["output_dir"], "bed", "all_cens.ct.bed"
        ),
    params:
        color=config["repeatmasker_sat_annot"]["ct_track_color"],
    log:
        "logs/create_ct_track.log",
    conda:
        "../env/tools.yaml"
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
            repeat=ANNOTATE_SAT_REPEATS.keys(),
        ),
        ct_track=rules.create_ct_track.output,
    output:
        os.path.join(
            config["repeatmasker_sat_annot"]["output_dir"],
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
            config["repeatmasker_sat_annot"]["output_dir"],
            "repeats",
            "all_cens_{chr}.annotation.fa.out",
        ),
    params:
        chr_pattern=lambda wc: str(wc.chr).replace("chr", "(chr|cen)") + "[_:\-\\tv]",
    shell:
        """
        grep -P "{params.chr_pattern}" {input.all_annotations} > {output.chr_annot}
        """


rule plot_satellite_annotations:
    input:
        script="workflow/scripts/repeatStructure_satellites.R",
        chr_annot=rules.split_rm_satellite_annotations.output,
    output:
        chr_plot=os.path.join(
            config["repeatmasker_sat_annot"]["output_dir"],
            "plots",
            "all_cens_{chr}.annotation.png",
        ),
    log:
        "logs/plot_{chr}_satellite_annotations.log",
    conda:
        "../env/r.yaml"
    shell:
        """
        Rscript {input.script} {input.chr_annot} {output.chr_plot} 2> {log}
        """


rule plot_repeatmasker_sat_only:
    input:
        expand(rules.plot_satellite_annotations.output, chr=CHROMOSOMES),
    default_target: True
