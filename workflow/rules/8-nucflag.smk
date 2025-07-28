include: "common.smk"
include: "5-ident_cen_ctgs.smk"
include: "7-fix_cens_w_repeatmasker.smk"


NUCFLAG_OUTDIR = join(OUTPUT_DIR, "8-nucflag")
NUCFLAG_LOGDIR = join(LOG_DIR, "8-nucflag")
NUCFLAG_BMKDIR = join(BMK_DIR, "8-nucflag")


if config.get("humas_annot"):

    include: "8-humas_annot.smk"


# Simplify RepeatMasker annotations
# Convert coords from relative to absolute coords.
# Add color specific for alpha-satellite.
rule format_repeatmasker_to_overlay_bed:
    input:
        rm=rules.fix_cens_rm_out.output,
    output:
        # [ctg_name, st, end, desc, action]
        join(
            NUCFLAG_OUTDIR,
            "{sm}_plot_rm.bed",
        ),
    log:
        join(NUCFLAG_LOGDIR, "format_repeatmasker_to_overlay_bed_{sm}.log"),
    params:
        color_alr_alpha="#8B008B",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        awk -v OFS="\\t" '{{
            name=$5; start=$6; end=$7; rType=$10; rClass=$11;

            # Find contig coordinates
            match(name, "^(.+):", ctg_name);
            match(name, ":(.+)-", ctg_start);
            match(name, ".*-(.+)$", ctg_end);
            new_start=start+ctg_start[1];
            new_end=end+ctg_start[1];

            # Split repeat class and replace specific repeat types.
            split(rClass, split_rClass, "/" );
            new_rClass=split_rClass[1];
            if (rClass == "Satellite/centr" || rClass == "Satellite") {{
                new_rClass=rType
            }}
            switch (new_rClass) {{
                case "SAR":
                    new_rClass="HSat1A";
                    break;
                case "HSAT":
                    new_rClass="HSat1B";
                    break;
                case "HSATII":
                    new_rClass="HSat2";
                    break;
                case "(CATTC)n":
                    new_rClass="HSat3";
                    break;
                case "(GAATG)n":
                    new_rClass="HSat3";
                    break;
                default:
                    break;
            }}

            # Set action for NucFlag
            action="plot"
            if (new_rClass == "ALR/Alpha") {{
                action="plot:{params.color_alr_alpha}"
            }}
            print ctg_name[1], new_start, new_end, new_rClass, action
        }}' {input.rm} > {output} 2> {log}
        """


# Convert stv_row bed to overlay bed with colors based on monomer len.
rule format_stv_to_overlay_bed:
    input:
        stv=lambda wc: humas_annot_sm_outputs(wc) if config.get("humas_annot") else [],
        annot_colors=(
            config["plot_hor_stv"]["stv_annot_colors"]
            if config.get("plot_hor_stv")
            else "config/stv_annot_colors.tsv"
        ),
    output:
        join(
            NUCFLAG_OUTDIR,
            "{sm}_plot_stv_row.bed",
        ),
    log:
        join(NUCFLAG_LOGDIR, "format_stv_to_overlay_bed_{sm}.log"),
    params:
        hor_mon_len=170,
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
            print ctg_name[1], st, end, num_mon, "plot:"stv_color
        }}' {input.annot_colors} {input.stv} > {output} 2> {log}
        """


# Create bedfile that only looks at regions annotated as ALR/Alpha and ignores everything else.
# Also add ignore bed if provided.
rule format_rm_nucflag_ignore_bed:
    input:
        bed=rules.format_repeatmasker_to_overlay_bed.output,
        ignore_bed=(
            config["nucflag"]["ignore_regions"]
            if config["nucflag"].get("ignore_regions")
            else []
        ),
    output:
        join(
            NUCFLAG_OUTDIR,
            "{sm}_correct_ALR_regions.rm.simple.bed",
        ),
    params:
        script=workflow.source_path("../scripts/simplify_rm_coords.py"),
        # Concatenate ignore bed.
        ignore_bed=lambda wc, input: (
            f"| cat - {input.ignore_bed}" if input.ignore_bed else ""
        ),
    conda:
        "../envs/py.yaml"
    log:
        join(NUCFLAG_LOGDIR, "simplify_rm_overlay_bed_{sm}.log"),
    shell:
        """
        python {params.script} -i {input.bed} {params.ignore_bed} > {output} 2> {log}
        """


# Create bedfile that only looks at live HOR and ignores everything else.
# Also add ignore bed if provided.
rule format_stv_nucflag_ignore_bed:
    input:
        stv=lambda wc: humas_annot_sm_outputs(wc) if config.get("humas_annot") else [],
        bed=ancient(rules.make_complete_cens_bed.output.cen_bed),
        ignore_bed=(
            config["nucflag"]["ignore_regions"]
            if config["nucflag"].get("ignore_regions")
            else []
        ),
    output:
        join(
            NUCFLAG_OUTDIR,
            "{sm}_ignore_stv_row.bed",
        ),
    params:
        bp_annot_gap_thr=1,
        # Concatenate ignore bed.
        ignore_bed=lambda wc, input: (
            f"| cat - {input.ignore_bed}" if input.ignore_bed else ""
        ),
    conda:
        "../envs/tools.yaml"
    log:
        join(NUCFLAG_LOGDIR, "format_stv_nucflag_ignore_bed_{sm}.log"),
    shell:
        """
        # Subtract all other regions from annotated HORs.
        # Include annotation gaps greater than {params.bp_annot_gap_thr} bp.
        {{ bedtools subtract \
        -a {input.bed} \
        -b <(cat {input.stv}) | \
        awk -v OFS="\\t" '{{
            len=$3-$2;
            if (len > {params.bp_annot_gap_thr}) {{
                print $0, "non-HOR", "ignore:absolute"
            }}
        }}' {params.ignore_bed};}} > {output} 2> {log}
        """


IGNORE_TYPE = config["nucflag"].get("ignore_type")
if IGNORE_TYPE == "asat":
    ignore_regions = str(rules.format_rm_nucflag_ignore_bed.output)
    overlay_beds = [str(rules.format_repeatmasker_to_overlay_bed.output)]
elif config.get("humas_annot") and IGNORE_TYPE == "live_asat":
    ignore_regions = [str(rules.format_stv_nucflag_ignore_bed.output)]
    overlay_beds = [str(rules.format_stv_to_overlay_bed.output)]
else:
    ignore_regions = str(rules.format_rm_nucflag_ignore_bed.output)
    overlay_beds = [str(rules.format_repeatmasker_to_overlay_bed.output)]

NUCFLAG_CFG = {
    "samples": [
        {
            "name": sm,
            "asm_fa": rules.rename_reort_asm.output.fa,
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
            "region_bed": ancient(rules.make_complete_cens_bed.output.cen_bed),
            # Ignore regions.
            "ignore_bed": ignore_regions,
            "overlay_beds": overlay_beds,
        }
        for sm in SAMPLE_NAMES
    ],
    "output_dir": NUCFLAG_OUTDIR,
    "logs_dir": NUCFLAG_LOGDIR,
    "benchmarks_dir": NUCFLAG_BMKDIR,
    **config["nucflag"],
}


module NucFlag:
    snakefile:
        "Snakemake-NucFlag/workflow/Snakefile"
    config:
        NUCFLAG_CFG


use rule * from NucFlag


rule nucflag_all:
    input:
        expand(rules.nucflag.input, sm=SAMPLE_NAMES),
    default_target: True
