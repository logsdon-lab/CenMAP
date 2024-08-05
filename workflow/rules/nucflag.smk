include: "common.smk"


rule format_repeatmasker_to_overlay_bed:
    input:
        rm=os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "{sm}",
            "{sm}_correct_ALR_regions.fa.reformatted.out",
        ),
    output:
        os.path.join(
            config["nucflag"]["output_dir"],
            "{sm}_correct_ALR_regions.rm.bed",
        ),
    log:
        "logs/nucflag/format_repeatmasker_to_overlay_bed_{sm}.log",
    params:
        color_alr_alpha="#8B008B",
    conda:
        "../env/tools.yaml"
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


rule simplify_rm_overlay_bed:
    input:
        script="workflow/scripts/simplify_rm_coords.py",
        bed=rules.format_repeatmasker_to_overlay_bed.output,
    output:
        os.path.join(
            config["nucflag"]["output_dir"],
            "{sm}_correct_ALR_regions.rm.simple.bed",
        ),
    params:
        color_alr_alpha="#8B008B",
        color_other="gray",
    conda:
        "../env/py.yaml"
    log:
        "logs/nucflag/simplify_rm_overlay_bed_{sm}.log",
    shell:
        """
        python {input.script} -i {input.bed} \
        --color_alr_alpha "{params.color_alr_alpha}" \
        --color_other "{params.color_other}" > {output} 2> {log}
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
                config["new_cens"]["output_dir"], "bed", f"{sm}_ALR_regions.bed"
            ),
            "overlay_beds": [
                # Original repeatmasker options
                str(rules.format_repeatmasker_to_overlay_bed.output),
                # Just ALR
                str(rules.simplify_rm_overlay_bed.output),
            ],
        }
        for sm in SAMPLE_NAMES
    ],
    **config["nucflag"],
}


module NucFlag:
    snakefile:
        github("logsdon-lab/Snakemake-NucFlag", path="workflow/Snakefile", tag="v0.1.0")
    config:
        NUCFLAG_CFG


use rule * from NucFlag


rule nucflag_only:
    input:
        expand(rules.format_repeatmasker_to_overlay_bed.output, sm=SAMPLE_NAMES),
        expand(rules.simplify_rm_overlay_bed.output, sm=SAMPLE_NAMES),
        expand(rules.nucflag.input, sm=SAMPLE_NAMES),
