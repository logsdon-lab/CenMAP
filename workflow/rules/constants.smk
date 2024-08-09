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

if config.get("plot_hor_stv"):
    with open(config["plot_hor_stv"]["sat_annot_colors"]) as fh:
        ANNOTATE_SAT_REPEATS = json.load(fh)
else:
    ANNOTATE_SAT_REPEATS = {}

REF_CENS_EDGE_LEN = round(
    (500_000 + config["extract_ref_hor_arrays"].get("added_bases", 0)) / 1000
)

# Thresholds
DEF_DNA_BRNN_FULL_ALR_THR = 30_000
DEF_CENSTATS_STATUS_EDGE_LEN_THR = 500_000
DEF_CENSTATS_STATUS_EDGE_PERC_ALR_THR = 0.95
DEF_CENSTATS_STATUS_MAX_ALR_LEN_THR = 250_000
