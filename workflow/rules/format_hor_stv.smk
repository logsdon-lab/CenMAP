include: "common.smk"
include: "humas_sd.smk"


checkpoint aggregate_format_all_stv_row:
    input:
        unpack(humas_sd_stv_outputs),
    output:
        os.path.join(
            config["plot_hor_stv"]["output_dir"], "bed", "{chr}_AS-HOR_stv_row.all.bed"
        ),
    log:
        "logs/plot_hor_stv/get_stv_row_from_{chr}_humas_hmmer_out.log",
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
        grep -P "{wildcards.chr}[_:]" || true ;}} > {output} 2> {log}
        """


rule filter_complete_correct_stv_row:
    input:
        stv_row_bed=os.path.join(
            config["plot_hor_stv"]["output_dir"], "bed", "{chr}_AS-HOR_stv_row.all.bed"
        ),
        complete_cens_bed=expand(
            os.path.join(
                config["ident_cen_ctgs"]["output_dir"],
                "bed",
                "final",
                "{sm}_complete_correct_cens.bed",
            ),
            sm=SAMPLE_NAMES,
        ),
    output:
        os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "bed",
            "{chr}_AS-HOR_stv_row.complete.bed",
        ),
    shell:
        """
        ( grep -f <(cut -f 1 {input.complete_cens_bed}) {input.stv_row_bed} || true ) > {output}
        """


# Get HOR monomer ort and merge monomers enforcing strandness.
# TODO: This should be doable in R but there are no correct interval libraries that meet all requirements:
# * are equivalent to bedtools without outright just wrapping bedtools (bedr, bedtoolsr, ...)
# * are correct (valr - removes interval edges when merging)
# * are simple/tidy (grange - wth)
# TODO: Remove R
rule get_stv_row_ort_bed:
    input:
        stv_row_bed=rules.aggregate_format_all_stv_row.output,
    output:
        # 4-col BED (chrom, start, end, strand)
        os.path.join(
            config["plot_hor_stv"]["output_dir"], "bed", "{chr}_AS-HOR_stv_row.ort.bed"
        ),
    params:
        dst_merge=100_000,
    conda:
        "../envs/tools.yaml"
    log:
        "logs/plot_cen_moddotplot/get_hor_mon_ort_{chr}.log",
    shell:
        """
        # Only get orientation of live HOR.
        bedtools merge -i <(sort -k1,1 -k2,2n {input} | awk -v OFS="\\t" '$4 ~ "L"') -s -d {params.dst_merge} -c 6 -o distinct > {output} 2> {log}
        """


rule format_hor_stv_only:
    input:
        expand(
            rules.aggregate_format_all_stv_row.output,
            chr=CHROMOSOMES,
        ),
        expand(
            rules.filter_complete_correct_stv_row.output,
            chr=CHROMOSOMES,
        ),
        expand(
            rules.get_stv_row_ort_bed.output,
            chr=CHROMOSOMES,
        ),