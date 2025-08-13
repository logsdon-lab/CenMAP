include: "common.smk"
include: "8-humas_annot.smk"
include: "9-get_complete_correct_cens.smk"


FMT_HOR_STV_OUTDIR = join(OUTPUT_DIR, "10-format_hor_stv")
FMT_HOR_STV_LOGDIR = join(LOG_DIR, "10-format_hor_stv")


checkpoint aggregate_format_all_stv_row:
    input:
        humas_annot_chr_outputs,
    output:
        join(FMT_HOR_STV_OUTDIR, "bed", "{chr}", "stv_all.bed"),
    log:
        join(FMT_HOR_STV_LOGDIR, "aggregate_format_all_stv_row_{chr}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ awk -v OFS="\\t" '{{
            # Find start in contig name.
            match($1, ":(.+)-", starts);
            if ($2 > $3) {{
                print "Invalid row at " NR, $0 > "{log}"
                next;
            }}
            # Update start and end
            $2=$2+starts[1];
            $3=$3+starts[1];
            $7=$7+starts[1];
            $8=$8+starts[1];
            print
        }}' {input} | \
        grep -P "{wildcards.chr}[_:-]" || true ;}} > {output} 2> {log}
        """


rule filter_complete_correct_stv_row:
    input:
        stv_row_bed=rules.aggregate_format_all_stv_row.output,
        complete_cens_bed=expand(
            rules.get_complete_correct_cens_bed.output,
            sm=SAMPLE_NAMES,
        ),
    output:
        join(
            FMT_HOR_STV_OUTDIR,
            "bed",
            "{chr}",
            "stv_complete.bed",
        ),
    shell:
        """
        ( grep -f <(cut -f 4 {input.complete_cens_bed} | awk '{{ print "^"$1":" }}') {input.stv_row_bed} || true ) > {output}
        """


rule format_hor_stv_all:
    input:
        expand(
            rules.aggregate_format_all_stv_row.output,
            chr=CHROMOSOMES,
        ),
        expand(
            rules.filter_complete_correct_stv_row.output,
            chr=CHROMOSOMES,
        ),
    default_target: True
