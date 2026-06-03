include: "common.smk"
include: "5-ident_cen_ctgs.smk"
include: "7-finalize_cens.smk"


NUCFLAG_OUTDIR = join(OUTPUT_DIR, "8-nucflag")
NUCFLAG_LOGDIR = join(LOG_DIR, "8-nucflag")
NUCFLAG_BMKDIR = join(BMK_DIR, "8-nucflag")
FILTER_BY_LIVE_ASAT = (
    config.get("humas_annot") and config["nucflag"].get("ignore_type") == "live_asat"
)

if config.get("humas_annot"):

    include: "8-humas_annot.smk"


# Simplify RepeatMasker annotations
# Convert coords from relative to absolute coords.
# Add color specific for alpha-satellite.
rule create_rm_overlay_bed:
    input:
        rm=(
            rules.create_fixed_rm_bed_by_sm.output
            if RUN_REPEATMASKER
            else rules.make_srf_putative_alr_regions.output
        ),
    output:
        # BED9
        join(
            NUCFLAG_OUTDIR,
            "{sm}_plot_rm.bed",
        ),
    log:
        join(NUCFLAG_LOGDIR, "create_rm_overlay_bed_{sm}.log"),
    params:
        script=(
            """<(printf '{{ match($1, "^(.+):", chroms); $1=chroms[1]; print }}')"""
            if RUN_REPEATMASKER
            else """<(printf '{{ print $1, $2, $3, $4, 0, ".", $2, $3, "#522758"}}')"""
        ),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        awk -v OFS="\\t" -f {params.script} {input.rm} > {output} 2> {log}
        """


# Convert stv_row bed to overlay bed with colors based on monomer len.
rule create_stv_overlay_bed:
    input:
        stv=lambda wc: humas_annot_sm_outputs(wc) if config.get("humas_annot") else [],
        # else branch should never be reached.
        annot_colors=(
            config["plot_hor_stv"]["stv_annot_colors"]
            if config.get("plot_hor_stv")
            else []
        ),
    output:
        join(
            NUCFLAG_OUTDIR,
            "{sm}_plot_stv_row.bed",
        ),
    log:
        join(NUCFLAG_LOGDIR, "create_stv_overlay_bed_{sm}.log"),
    params:
        hor_mon_len=170,
        stv=lambda wc, input: input.stv if input.stv else '<(echo "")',
    conda:
        "../envs/tools.yaml"
    shell:
        """
        awk -v OFS="\\t" '{{
            # Load stv colors from first file.
            if (NR == FNR) {{ stv_colors[$1]=$2; next;}}
            name=$1;
            match(name, "^(.+):", ctg_name);
            st=$2; end=$3;
            hor=$4;
            hor_len=$3-$2;
            num_mon=int(hor_len / {params.hor_mon_len})
            stv_color=stv_colors[num_mon]
            if (stv_color == "") {{ stv_color="gray"; }}
            print ctg_name[1], st, end, num_mon, 0, ".", st, end, stv_color
        }}' {input.annot_colors} {params.stv} > {output} 2> {log}
        """


# Create bedfile that only looks at live asat HOR regions or asat regions.
rule create_nucflag_ignore_bed:
    input:
        bed9_annot=lambda wc: (
            humas_annot_sm_outputs(wc)
            if FILTER_BY_LIVE_ASAT
            else rules.create_rm_overlay_bed.output
        ),
        bed_cen=rules.make_complete_cens_bed.output.cen_bed,
    output:
        join(
            NUCFLAG_OUTDIR,
            "{sm}_ignore_regions.bed",
        ),
    params:
        # Size of LINE
        bp_merge=8000,
        filter_str="L" if FILTER_BY_LIVE_ASAT else "ALR",
        bed9_annot=lambda wc, input: (
            input.bed9_annot if input.bed9_annot else '<(echo "")'
        ),
    conda:
        "../envs/tools.yaml"
    log:
        join(NUCFLAG_LOGDIR, "create_nucflag_group_bed_{sm}.log"),
    shell:
        """
        # Filter to only live or asat. Sort
        # Merge by small distance
        # Then subtract by everything else
        {{ awk -v OFS="\\t" '{{if ($4 ~ "{params.filter_str}") {{print $1, $2, $3}}}}' {params.bed9_annot} | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i - -d {params.bp_merge} | \
        bedtools subtract \
            -a <(sort -k1,1 -k2,2n {input.bed_cen}) \
            -b - | \
        cut -f1-3 ;}} > {output} 2> {log}
        """


NUCFLAG_CFG = {
    "samples": [
        {
            "name": sm,
            "asm_fa": rules.create_final_asm.output.fa,
            # Switch between fofn dir or read dir + ext.
            **(
                {
                    "read_fofn": join(
                        config["nucflag"]["input_hifi_reads_fofn_dir"], f"{sm}.fofn"
                    ),
                }
                if config["nucflag"].get("input_hifi_reads_fofn_dir")
                else {
                    "read_dir": join(config["nucflag"]["input_hifi_reads_dir"], sm),
                    "read_rgx": config["nucflag"]["reads_rgx"],
                }
            ),
            "config": config["nucflag"]["config_nucflag"],
            "region_bed": rules.make_complete_cens_bed.output.cen_bed,
            # Ignore regions.
            "ignore_bed": [str(rules.create_nucflag_ignore_bed.output)],
            "overlay_beds": [
                str(
                    rules.create_stv_overlay_bed.output
                    if FILTER_BY_LIVE_ASAT
                    else rules.create_rm_overlay_bed.output
                )
            ],
        }
        for sm in SAMPLE_NAMES
    ],
    "output_dir": NUCFLAG_OUTDIR,
    "logs_dir": NUCFLAG_LOGDIR,
    "benchmarks_dir": NUCFLAG_BMKDIR,
    "output_plots": True,
    **config["nucflag"],
}


module NucFlag:
    snakefile:
        "Snakemake-NucFlag/workflow/Snakefile"
    config:
        NUCFLAG_CFG


use rule * from NucFlag


# NucFlag's default output status merges by non-overlapping regions.
# We want just by input region. Since we ignore non-asat, we do it here.
rule create_region_status:
    input:
        bed_cen=rules.make_complete_cens_bed.output.cen_bed,
        bed_calls=rules.check_asm_nucflag.output.misassemblies,
    output:
        join(
            NUCFLAG_OUTDIR,
            "{sm}_status_asat.bed",
        ),
    conda:
        "Snakemake-NucFlag/workflow/env/nucflag.yaml"
    log:
        join(NUCFLAG_LOGDIR, "create_region_status_{sm}.log"),
    shell:
        """
        nucflag status -i {input.bed_calls} -b {input.bed_cen} > {output} 2> {log}
        """


# Filter if:
# * Contains misjoin, deletion, insertion, other_repeat, false_duplication, collapse, scaffold
# * QV less than some value
# If no conditions, no operation done.
rule filter_nucflag_status_bed:
    input:
        rules.create_region_status.output,
    output:
        join(
            NUCFLAG_OUTDIR,
            "{sm}_status_asat_filtered.bed",
        ),
    params:
        script=workflow.source_path("../scripts/filter_nucflag_status.py"),
        lt_qv=(
            f"-q {config["nucflag"]["filter_if_lt_qv"]}"
            if config["nucflag"]["filter_if_lt_qv"]
            else ""
        ),
        contains_type=f"-t {" ".join(config["nucflag"]["filter_if_contains"])}",
    shell:
        """
        python {params.script} {input} {params.lt_qv} {params.contains_type} > {output}
        """


rule nucflag_all:
    input:
        expand(rules.nucflag.input, sm=SAMPLE_NAMES),
        expand(rules.filter_nucflag_status_bed.input, sm=SAMPLE_NAMES),
    default_target: True
