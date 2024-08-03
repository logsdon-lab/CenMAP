include: "common.smk"


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
            "ignore_bed": config["nucflag"]["ignore_regions"],
            "region_bed": os.path.join(
                config["new_cens"]["output_dir"], "bed", f"{sm}_ALR_regions.bed"
            ),
        }
        for sm in SAMPLE_NAMES
    ],
    **config["nucflag"],
}


module NucFlag:
    snakefile:
        github(
            "logsdon-lab/Snakemake-NucFlag", path="workflow/Snakefile", branch="main"
        )
    config:
        NUCFLAG_CFG


use rule * from NucFlag


rule nucflag_only:
    input:
        expand(rules.nucflag.input, sm=SAMPLE_NAMES),
