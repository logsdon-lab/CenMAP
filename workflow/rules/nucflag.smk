include: "common.smk"
include: "humas_sd.smk"


rule format_repeatmasker_to_overlay_bed:
    input:
        rm=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "all",
            "reoriented_{sm}_cens.fa.out",
        ),
    output:
        os.path.join(
            config["nucflag"]["output_dir"],
            "{sm}_plot_rm.bed",
        ),
    log:
        "logs/nucflag/format_repeatmasker_to_overlay_bed_{sm}.log",
    params:
        color_alr_alpha="#8B008B",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        awk -v OFS="\\t" '{{
            name=$5; start=$6; end=$7; rType=$10; rClass=$11;

            # Find contig coordinates
            match(name, "^(.+):", abbr_name);
            match(name, ":(.+)-", ctg_start);
            match(name, ".*-(.+)$", ctg_end);
            new_name=abbr_name[1];
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
                    new_rClass="HSat2";
                    break;
                case "(GAATG)n":
                    new_rClass="HSat2";
                    break;
                default:
                    break;
            }}

            # Set action for NucFlag
            action="plot"
            if (new_rClass == "ALR/Alpha") {{
                action="plot:{params.color_alr_alpha}"
            }}
            print new_name, new_start, new_end, new_rClass, action
        }}' {input.rm} > {output} 2> {log}
        """


rule format_stv_to_overlay_bed:
    input:
        stv=humas_sd_stv_sm_outputs,
    output:
        os.path.join(
            config["nucflag"]["output_dir"],
            "{sm}_plot_stv_row.bed",
        ),
    log:
        "logs/nucflag/format_stv_to_overlay_bed_{sm}.log",
    params:
        hor_mon_len=170,
    conda:
        "../envs/tools.yaml"
    shell:
        """
        awk -v OFS="\\t" '{{
            name=$1;
            st=$2; end=$3;
            hor=$4;
            hor_len=$3-$2;
            num_mon=int(hor_len / {params.hor_mon_len})
            switch (num_mon) {{
                case 1: color="\\#A8275C"; break;
                case 10: color="\\#9AC78A"; break;
                case 11: color="\\#CC8FC1"; break;
                case 12: color="\\#3997C6"; break;
                case 13: color="\\#8882C4"; break;
                case 14: color="\\#8ABDD6"; break;
                case 15: color="\\#096858"; break;
                case 16: color="\\#45B4CE"; break;
                case 17: color="\\#AFA7D8"; break;
                case 18: color="\\#A874B5"; break;
                case 19: color="\\#3F66A0"; break;
                case 2: color="\\#D66C54"; break;
                case 20: color="\\#BFDD97"; break;
                case 21: color="\\#AF5D87"; break;
                case 22: color="\\#E5E57A"; break;
                case 24: color="\\#ED975D"; break;
                case 26: color="\\#F9E193"; break;
                case 3: color="\\#93430C"; break;
                case 30: color="\\#E5D1A1"; break;
                case 32: color="\\#A1B5E5"; break;
                case 34: color="\\#9F68A5"; break;
                case 35: color="\\#81B25B"; break;
                case 4: color="\\#F4DC78"; break;
                case 5: color="\\#7EC0B3"; break;
                case 6: color="\\#B23F73"; break;
                case 7: color="\\#8CC49F"; break;
                case 8: color="\\#893F89"; break;
                case 9: color="\\#6565AA"; break;
                default: color="gray"; break;
            }}
            print name, st, end, num_mon, "plot:"color
        }}' {input} > {output} 2> {log}
        """


# Create bedfile that only looks at live HOR and ignores everything else.
rule format_stv_nucflag_ignore_bed:
    input:
        fai=os.path.join(
            config["concat_asm"]["output_dir"],
            "{sm}",
            "{sm}_regions.renamed.reort.fa.fai",
        ),
        stv=humas_sd_stv_sm_outputs,
    output:
        os.path.join(
            config["nucflag"]["output_dir"],
            "{sm}_ignore_stv_row.bed",
        ),
    params:
        bp_annot_gap_thr=1,
    conda:
        "../envs/tools.yaml"
    log:
        "logs/nucflag/format_stv_nucflag_ignore_bed_{sm}.log",
    shell:
        """
        # Subtract all other regions from annotated HORs.
        # Include annotation gaps greater than {params.bp_annot_gap_thr} bp.
        {{ bedtools subtract \
        -a <(awk -v OFS="\\t" '{{ print $1, "1", $2 }}' {input.fai}) \
        -b {input.stv} | \
        awk -v OFS="\\t" '{{
            len=$3-$2;
            if (len > {params.bp_annot_gap_thr}) {{
                print $0, "non-HOR", "ignore:absolute"
            }}
        }}';}} > {output} 2> {log}
        """


NUCFLAG_CFG = {
    "samples": [
        {
            "name": sm,
            "asm_fa": os.path.join(
                config["concat_asm"]["output_dir"], sm, f"{sm}_regions.renamed.reort.fa"
            ),
            # Switch between fofn dir or read dir + ext.
            **(
                {
                    "read_fofn": os.path.join(
                        config["nucflag"]["hifi_reads_fofn_dir"], f"{sm}.fofn"
                    ),
                }
                if config["nucflag"].get("hifi_reads_fofn_dir")
                else {
                    "read_dir": os.path.join(config["nucflag"]["hifi_reads_dir"], sm),
                    "read_ext": config["nucflag"]["reads_ext"],
                }
            ),
            "config": config["nucflag"]["config_nucflag"],
            "region_bed": os.path.join(
                config["ident_cen_ctgs"]["output_dir"],
                "bed",
                "interm",
                "{sm}_complete_cens.bed",
            ),
            # # Ignore regions.
            "ignore_bed": str(rules.format_stv_nucflag_ignore_bed.output),
            "overlay_beds": [
                # Original repeatmasker options
                str(rules.format_repeatmasker_to_overlay_bed.output),
                str(rules.format_stv_to_overlay_bed.output),
            ],
        }
        for sm in SAMPLE_NAMES
    ],
    **config["nucflag"],
}


module NucFlag:
    snakefile:
        github("logsdon-lab/Snakemake-NucFlag", path="workflow/Snakefile", tag="v0.2.2")
    config:
        NUCFLAG_CFG


use rule * from NucFlag


rule nucflag_only:
    input:
        expand(rules.format_repeatmasker_to_overlay_bed.output, sm=SAMPLE_NAMES),
        expand(rules.format_stv_to_overlay_bed.output, sm=SAMPLE_NAMES),
        expand(rules.nucflag.input, sm=SAMPLE_NAMES),
