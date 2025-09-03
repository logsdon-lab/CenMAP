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
REF_NAME = "T2T-CHM13"
REF_URL = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
REF_FA = join(OUTPUT_DIR, "1-download_ref", f"{REF_NAME}.fa.gz")

SAMPLE_NAMES = config["samples"]
CHROMOSOMES = config["chromosomes"]
RGX_CHR = re.compile(r"(chr[0-9XY]+)")

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
