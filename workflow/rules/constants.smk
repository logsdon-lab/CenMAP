# Globals shared throughout workflow.
REF_FA = config["align_asm_to_ref"]["reference"]
REF_NAME = os.path.splitext(os.path.basename(REF_FA).strip(".gz"))[0]
SAMPLE_NAMES = config["samples"]
HAPLOTYPE = ("haplotype1", "haplotype2")
MER_ORDER = ("small", "large")
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
RGX_CHR = re.compile("(chr[0-9XY]+)")

with open(config["repeatmasker_sat_annot"]["config_pattern_colors"]) as fh:
    ANNOTATE_SAT_REPEATS = json.load(fh)

REF_CENS_EDGE_LEN = round(
    (500_000 + config["extract_ref_hor_arrays"].get("added_bases", 0)) / 1000
)
