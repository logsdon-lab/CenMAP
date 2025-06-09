# Globals shared throughout workflow.

# Reference genome.
REF_NAME = "T2T-CHM13"
REF_URL = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
REF_FA = f"data/reference/{REF_NAME}.fa.gz"

# Container for singularity
CONTAINER = config["container"]

# Dirs
OUTPUT_DIR = "results"
LOG_DIR = "logs"
BMK_DIR = "benchmarks"
HUMAS_CENS_SPLIT_DIR = join(OUTPUT_DIR, "8-humas_annot", "seq")
FINAL_OUTPUT_DIR = join(OUTPUT_DIR, "final")

SAMPLE_NAMES = config["samples"]
CHROMOSOMES = config["chromosomes"]
RGX_CHR = re.compile(r"(chr[0-9XY]+)")

REF_CENS_EDGE_LEN = round(
    (500_000 + config["extract_ref_hor_arrays"]["added_bases"]) / 1000
)

# Check if command is to containerize workflow.
ARGUMENTS = set(sys.argv)
IS_CONTAINERIZE_CMD = "--containerize" in ARGUMENTS
IS_SINGULARITY = (
    DeploymentMethod.APPTAINER in workflow.deployment_settings.deployment_method
)
IS_CONDA = DeploymentMethod.CONDA in workflow.deployment_settings.deployment_method


# Thresholds
DEF_CENSTATS_STATUS_EDGE_LEN_THR = 500_000
DEF_CENSTATS_STATUS_EDGE_PERC_ALR_THR = 0.95
DEF_CENSTATS_STATUS_MAX_ALR_LEN_THR = 250_000


# Set shared constraints.
wildcard_constraints:
    chr="|".join(CHROMOSOMES),
    sm="|".join(SAMPLE_NAMES),
