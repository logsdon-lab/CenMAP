import sys
import argparse
import pandas as pd


RM_COLS = [
    "idx",
    "div",
    "deldiv",
    "insdiv",
    "contig",
    "start",
    "end",
    "left",
    "C",
    "type",
    "rClass",
    "right",
    "x",
    "y",
    "z",
]
DEFAULT_EDGE_LEN = 500_000


def jaccard_index(a: set[str], b: set[str]) -> float:
    """
    Jaccard similarity index.
    * https://www.statisticshowto.com/jaccard-index/
    """
    return (len(a.intersection(b)) / len(a.union(b))) * 100.0


def is_reversed_cens(
    ref: pd.DataFrame, comp: pd.DataFrame, *, edge_len: int = DEFAULT_EDGE_LEN
) -> bool:
    """
    Check if centromere RepeatMasker output is in the reverse orientation with respect to a reference centromere RepeatMasker output.
    Determined via Jaccard similarity index computed from the unique repeat types at the edges of the central AS-HOR array.

    ### Params
    * `ref`: Reference centromere RepeatMasker output.
    * `comp`: Comparison centromere RepeatMasker output.
    * `edge_len`: Edge length to evaluate.

    ### Returns
    * If `comp` centromeres should be reversed.
    """
    ref_ledge = ref.loc[ref["end"] < edge_len]
    ref_redge = ref.loc[ref["end"] > ref.iloc[-1]["end"] - edge_len]

    comp_ledge = comp.loc[comp["end"] < edge_len]
    comp_redge = comp.loc[comp["end"] > comp.iloc[-1]["end"] - edge_len]

    ref_redge_types = set(ref_redge["type"].unique())
    ref_ledge_types = set(ref_ledge["type"].unique())

    comp_redge_types = set(comp_redge["type"].unique())
    comp_ledge_types = set(comp_ledge["type"].unique())

    rl_cr = jaccard_index(ref_ledge_types, comp_redge_types)
    rr_cl = jaccard_index(ref_redge_types, comp_ledge_types)
    rr_cr = jaccard_index(ref_redge_types, comp_redge_types)
    rl_cl = jaccard_index(ref_ledge_types, comp_ledge_types)

    highest_similarity_score = max((rl_cr, rr_cl, rr_cr, rl_cl))
    return highest_similarity_score in (rl_cr, rr_cl)


def main():
    ap = argparse.ArgumentParser(
        description="Determines if centromeres are incorrectly oriented with respect to a reference."
    )
    ap.add_argument(
        "-i",
        "--input",
        help="Input RepeatMasker output. Should contain contig reference. Expects no header.",
        default=sys.stdin,
        type=argparse.FileType("rt"),
    )
    ap.add_argument(
        "-o",
        "--output",
        help="List of contigs that are reversed.",
        default=sys.stdout,
        type=argparse.FileType("wt"),
    )
    ap.add_argument(
        "--reference_prefix", help="Reference contig name prefix.", default="chm13"
    )
    ap.add_argument(
        "--edge_len",
        help="Edge length to along contigs to evaluate.",
        default=DEFAULT_EDGE_LEN,
    )
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t", header=None, names=RM_COLS)
    ref_ctg_name = next(
        (ctg for ctg in df["contig"].unique() if ctg.startswith(args.reference_prefix))
    )
    df_ref = df.loc[df["contig"] == ref_ctg_name]
    for ctg in df["contig"].unique():
        df_ctg = df.loc[df["contig"] == ctg]
        if is_reversed_cens(ref=df_ref, comp=df_ctg):
            args.output.write(f"{ctg}\n")


if __name__ == "__main__":
    raise SystemExit(main())
