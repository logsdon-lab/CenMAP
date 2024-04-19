# Globals shared throughout workflow.
REF_FA = config["align_asm_to_ref"]["reference"]
REF_NAME = os.path.splitext(os.path.basename(REF_FA))[0]
SAMPLE_NAMES = config["samples"]
ORIENTATION = ("fwd", "rev")
HAPLOTYPE = ("haplotype1", "haplotype2")
CHROMOSOMES = config.get(
    "chromosomes",
    (
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY",
    ),
)
SAMPLE_FLOWCELL_IDS = get_hifi_read_wildcards()

with open(config["repeatmasker_sat_annot"]["config_pattern_colors"]) as fh:
    ANNOTATE_SAT_REPEATS = json.load(fh)
