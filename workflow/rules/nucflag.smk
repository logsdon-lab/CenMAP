NUCFLAG_CFG = {
    "samples": [
        {
            "name": sm,
            "asm_fa": rules.concat_asm.output,
            "read_dir": os.path.join(config["nucflag"]["hifi_reads_dir"], sm),
            "read_ext": config["nucflag"]["reads_ext"],
            "region_bed": rules.make_new_cens_bed_file.output.alr_bed,
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
