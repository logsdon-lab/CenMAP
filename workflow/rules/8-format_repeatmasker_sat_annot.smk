include: "common.smk"
include: "7-fix_cens_w_repeatmasker.smk"


FMT_RM_SAT_OUTDIR = join(OUTPUT_DIR, "8-format_repeatmasker_sat_annot")
FMT_RM_SAT_LOGDIR = join(LOG_DIR, "8-format_repeatmasker_sat_annot")
FMT_RM_SAT_BMKDIR = join(BMK_DIR, "8-format_repeatmasker_sat_annot")


rule merge_complete_and_correct_rm_out:
    input:
        expand(
            rules.fix_cens_rm_out.output,
            sm=SAMPLE_NAMES,
        ),
    output:
        join(
            FMT_RM_SAT_OUTDIR,
            "repeats",
            "all",
            "reoriented_all_samples_and_ref_cens.fa.out",
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
        join(
            FMT_RM_SAT_OUTDIR,
            "bed",
            "all_cens_{repeat}.satellites.bed",
        ),
    conda:
        "../envs/tools.yaml"
    params:
        pattern=lambda wc: ANNOTATE_SAT_REPEATS.get(str(wc.repeat), {}).get("pattern"),
        color=lambda wc: ANNOTATE_SAT_REPEATS.get(str(wc.repeat), {}).get("color"),
    log:
        join(FMT_RM_SAT_LOGDIR, "create_{repeat}_annotated_satellites.log"),
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
        corrected_rm_out=rules.merge_complete_and_correct_rm_out.output,
    output:
        join(FMT_RM_SAT_OUTDIR, "bed", "all_cens.ct.bed"),
    params:
        color=ANNOTATE_SAT_REPEATS.get("ct", {}).get("color"),
    log:
        join(FMT_RM_SAT_LOGDIR, "create_ct_track.log"),
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
        join(
            FMT_RM_SAT_OUTDIR,
            "bed",
            "all_cens.annotation.bed",
        ),
    shell:
        """
        cat {input.ct_track} {input.satellites} > {output}
        """


rule format_repeatmasker_sat_all:
    input:
        rules.aggregate_rm_satellite_annotations.output,
    default_target: True
