# Globals shared throughout workflow.

# Container for singularity
CONTAINER = config["container"]

# Dirs
OUTPUT_DIR = config["output_dir"]
LOG_DIR = config["log_dir"]
BMK_DIR = config["benchmark_dir"]
HUMAS_CENS_SPLIT_DIR = join(OUTPUT_DIR, "8-humas_annot", "seq")
FINAL_OUTPUT_DIR = join(OUTPUT_DIR, "final")

# Reference genome.
CHROMOSOMES = config.get("chromosomes", [])
if config["ident_cen_ctgs"].get("reference"):
    REF_NAME, _, _ = os.path.splitext(
        os.path.basename(config["ident_cen_ctgs"]["reference"])
    )[0].partition(".")
    REF_URL = None
    REF_FA = config["ident_cen_ctgs"]["reference"]
# No reference at all. Don't stratify output by chromosomes.
elif not CHROMOSOMES:
    REF_NAME = None
    REF_URL = None
    REF_FA = []
else:
    REF_NAME = "T2T-CHM13"
    REF_URL = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
    REF_FA = join(OUTPUT_DIR, "1-download_ref", f"{REF_NAME}.fa.gz")

SAMPLE_NAMES = config["samples"]
IS_HUMAN_ANNOT = config.get("humas_annot", {}).get("mode") in ("sd", "hmmer")
RUN_REPEATMASKER = "repeatmasker" in config

RGX_CHR = re.compile(r"(chr[0-9XY]+)")
RGX_SM_CTG = re.compile(r"^(.+)_(.+)$")
# This monstrosity matches the expected patterns.
# https://regex101.com/r/u4A2Xj/1
RGX_SM_CHR_CTG = re.compile(
    r"^(.*?)_((?:(?:(?:chr|rc-chr)[0-9XY]+)-)*(?:(?:chr|rc-chr)[0-9XY]+))_(.+)$"
)

# Check if command is to containerize workflow.
ARGUMENTS = set(sys.argv)
IS_CONTAINERIZE_CMD = "--containerize" in ARGUMENTS
IS_SINGULARITY = (
    DeploymentMethod.APPTAINER in workflow.deployment_settings.deployment_method
)
IS_CONDA = DeploymentMethod.CONDA in workflow.deployment_settings.deployment_method


# Set shared constraints.
wildcard_constraints:
    chr="|".join(CHROMOSOMES),
    sm="|".join(SAMPLE_NAMES),
