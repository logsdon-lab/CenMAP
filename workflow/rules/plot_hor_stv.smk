include: "common.smk"
include: "humas_hmmer.smk"


rule get_live_hor:
    input:
        script="workflow/scripts/stv_fix/scripts/live_HORs_filter.py",
        # fname is wildcard
        humas_hmmer_out=rules.cens_filter_hmm_res_overlaps_as_hor.output,
    output:
        renamed_bed=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "bed",
            "{chr}",
            "{fname}_renamed.bed",
        ),
        live_hor_bed=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "bed",
            "{chr}",
            "{fname}_liveHORs.bed",
        ),
    params:
        # Fix inversion for chromosome 1 and 19.
        inversion_fix_cmd=lambda wc: (
            "| sed 's+S1C1/5/19H1L.6/4+S1C1/5/19H1L.6+g'"
            if wc.chr == "chr1" or wc.chr == "chr19"
            else ""
        ),
    log:
        "logs/plot_hor_stv/get_live_hor_{fname}_{chr}.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        awk -v OFS="\\t" '{{
            print "{wildcards.chr}", $2, $3, $4, $5, $6, $7, $8, $9, $1
        }}' {input.humas_hmmer_out} > {output.renamed_bed} 2> {log}
        {{ python3 {input.script} {output.renamed_bed} {params.inversion_fix_cmd} ;}} > {output.live_hor_bed} 2>> {log}
        """


rule filter_as_hor_stv_bed:
    input:
        script="workflow/scripts/stv_fix/scripts/mon2stv.py",
        live_hor_bed=rules.get_live_hor.output.live_hor_bed,
        renamed_bed=rules.get_live_hor.output.renamed_bed,
    output:
        stv_row_bed=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "bed",
            "{chr}",
            "{fname}_stv_row.bed",
        ),
        as_hor_stv_row_bed=os.path.join(
            config["plot_hor_stv"]["output_dir"],
            "bed",
            "{chr}",
            "AS-HOR_{fname}_stv_row.bed",
        ),
    log:
        "logs/plot_hor_stv/filter_as_hor_stv_bed_{fname}_{chr}.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        python3 {input.script} {input.live_hor_bed} > {output.stv_row_bed} 2>> {log}
        awk -v OFS="\\t" 'FNR==NR{{a[NR]=$10;next}}{{$1=a[FNR]}}1' {input.renamed_bed} {output.stv_row_bed} > {output.as_hor_stv_row_bed} 2>> {log}
        """


def as_hor_bedfiles(wc):
    _ = checkpoints.split_cens_for_humas_hmmer.get(**wc).output
    fnames, chrs = extract_fa_fnames_and_chr(
        config["humas_hmmer"]["input_dir"], filter_chr=str(wc.chr)
    )
    _ = checkpoints.run_humas_hmmer_for_anvil.get(**wc).output
    return dict(
        stv_bed_filtered=expand(
            rules.filter_as_hor_stv_bed.output, zip, fname=fnames, chr=chrs
        )
    )


checkpoint aggregate_format_all_stv_row:
    input:
        unpack(as_hor_bedfiles),
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
        awk -v OFS="\\t" '{{
            # Find start in contig name.
            match($1, ":(.+)-", starts);
            # Update start and end
            $2=$2+starts[1];
            $3=$3+starts[1];
            $7=$7+starts[1];
            $8=$8+starts[1];
            print
        }}' {input} > {output}
        """


rule plot_stv_with_order:
    input:
        script="workflow/scripts/plot_cens_stvHOR.R",
        chm1_stv=config["plot_hor_stv"]["chm1_stv"],
        chm13_stv=config["plot_hor_stv"]["chm13_stv"],
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
            --input_chm13 {input.chm13_stv} \
            --input_chm1 {input.chm1_stv} \
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
