include: "common.smk"
include: "humas_hmmer.smk"


checkpoint aggregate_format_all_stv_row:
    input:
        unpack(humas_hmmer_stv_outputs),
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
        grep -P "{wildcards.chr}[_:]" ;}} > {output} 2> {log}
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
        bedtools merge -i <(sort -k1,1 -k2,2n {input}) -s -d {params.dst_merge} -c 6 -o distinct > {output} 2> {log}
        """


# No ort added.
rule plot_stv_with_order:
    input:
        script="workflow/scripts/plot_cens_stvHOR.R",
        all_stv=rules.aggregate_format_all_stv_row.output,
    output:
        hor_array_plot=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "plots",
            "hor",
            "{chr}.png",
        ),
    params:
        mer_order=lambda wc: MONOMER_ORDER[wc.chr],
    log:
        "logs/plot_hor_stv/plot_{chr}_stv.log",
    conda:
        "../envs/r.yaml"
    shell:
        """
        if ! [ -s {input.all_stv} ]; then
            touch {output.hor_array_plot}
        else
            Rscript {input.script} \
            --input {input.all_stv} \
            --output {output.hor_array_plot} \
            --chr {wildcards.chr} \
            --mer_order {params.mer_order} 2> {log}
        fi
        """


rule plot_hor_stv_only:
    input:
        expand(
            rules.plot_stv_with_order.output,
            chr=CHROMOSOMES,
        ),
