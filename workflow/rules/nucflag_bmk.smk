include: "common.smk"
include: "download_data.smk"


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


rule simplify_rm_overlay_bed:
    input:
        script="workflow/scripts/simplify_rm_coords.py",
        bed=rules.format_repeatmasker_to_overlay_bed.output,
    output:
        os.path.join(
            config["nucflag"]["output_dir"],
            "{sm}_correct_ALR_regions.rm.simple.bed",
        ),
    conda:
        "../envs/py.yaml"
    log:
        "logs/nucflag/simplify_rm_overlay_bed_{sm}.log",
    shell:
        """
        python {input.script} -i {input.bed} > {output} 2> {log}
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
            # Ignore regions.
            "ignore_bed": str(rules.simplify_rm_overlay_bed.output),
            "overlay_beds": [
                # Original repeatmasker options
                str(rules.format_repeatmasker_to_overlay_bed.output),
            ],
        }
        for sm in SAMPLE_NAMES
    ],
    **config["nucflag"],
}

import copy

WINNOWMAP_CFG = copy.deepcopy(NUCFLAG_CFG)
MINIMAP2_CFG = copy.deepcopy(NUCFLAG_CFG)
PBMM2_CFG = copy.deepcopy(NUCFLAG_CFG)


def update_cfg_aligner(cfg: dict, aligner: str) -> dict:
    cfg["aligner"] = aligner
    cfg["output_dir"] = f"results/nucflag_{aligner}"
    cfg["tmp_dir"] = f"temp_{aligner}"
    cfg["logs_dir"] = f"logs/nucflag_{aligner}"
    cfg["benchmarks_dir"] = f"benchmarks/nucflag_{aligner}"
    return cfg


WINNOWMAP_CFG = update_cfg_aligner(WINNOWMAP_CFG, "winnowmap")
MINIMAP2_CFG = update_cfg_aligner(MINIMAP2_CFG, "minimap2")
PBMM2_CFG = update_cfg_aligner(PBMM2_CFG, "pbmm2")


module NucFlag_Winnowmap:
    snakefile:
        github(
            "logsdon-lab/Snakemake-NucFlag",
            path="workflow/Snakefile",
            branch="feature/choose-aligner",
        )
    config:
        WINNOWMAP_CFG


module NucFlag_Minimap2:
    snakefile:
        github(
            "logsdon-lab/Snakemake-NucFlag",
            path="workflow/Snakefile",
            branch="feature/choose-aligner",
        )
    config:
        MINIMAP2_CFG


module NucFlag_Pbmm2:
    snakefile:
        github(
            "logsdon-lab/Snakemake-NucFlag",
            path="workflow/Snakefile",
            branch="feature/choose-aligner",
        )
    config:
        PBMM2_CFG


use rule * from NucFlag_Winnowmap as wm_*


use rule * from NucFlag_Minimap2 as mm2_*


use rule * from NucFlag_Pbmm2 as pbmm2_*


rule summarize_benchmark:
    input:
        expand(rules.wm_nucflag.input, sm=SAMPLE_NAMES),
        expand(rules.mm2_nucflag.input, sm=SAMPLE_NAMES),
        expand(rules.pbmm2_nucflag.input, sm=SAMPLE_NAMES),
    output:
        plt_cens_status=os.path.join(
            "results", "nucflag_aligner_comparison", "cens_status.png"
        ),
        plt_mem_usage=os.path.join(
            "results", "nucflag_aligner_comparison", "memory_usage.png"
        ),
        plt_cens_venn=os.path.join(
            "results", "nucflag_aligner_comparison", "venn_shared_correct_cens.png"
        ),
        plt_wall_time=os.path.join(
            "results", "nucflag_aligner_comparison", "wall_time.png"
        ),
    conda:
        "../envs/jupyter.yaml"
    notebook:
        "../notebooks/summarize_bmks.ipynb"


rule nucflag_bmk_only:
    input:
        expand(rules.download_bmk_data_all.input, sm=SAMPLE_NAMES),
        expand(rules.format_repeatmasker_to_overlay_bed.output, sm=SAMPLE_NAMES),
        expand(rules.simplify_rm_overlay_bed.output, sm=SAMPLE_NAMES),
        expand(rules.wm_nucflag.input, sm=SAMPLE_NAMES),
        expand(rules.mm2_nucflag.input, sm=SAMPLE_NAMES),
        expand(rules.pbmm2_nucflag.input, sm=SAMPLE_NAMES),
        rules.summarize_benchmark.output,
    default_target: True
