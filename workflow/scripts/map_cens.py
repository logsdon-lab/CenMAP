import sys
import argparse
import polars as pl

ALN_HEADER = [
    "#reference_name",
    "reference_start",
    "reference_end",
    "reference_length",
    "strand",
    "query_name",
    "query_start",
    "query_end",
    "query_length",
    "perID_by_matches",
    "perID_by_events",
    "perID_by_all",
    "matches",
    "mismatches",
    "deletion_events",
    "insertion_events",
    "deletions",
    "insertions",
]


def main():
    ap = argparse.ArgumentParser("Map cens from alignments.")
    ap.add_argument("-i", "--input", help="Input alignment bed.", type=str)
    ap.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))
    ap.add_argument("-t", "--len_thr", default=1_000_000, type=int)

    args = ap.parse_args()

    df = pl.read_csv(
        args.input, separator="\t", new_columns=ALN_HEADER, has_header=False
    )

    # Get concensus mapping by choosing reference chr with most rows mapping to it.
    df_concensus_mapping = (
        df.group_by(["#reference_name", "query_name", "strand"])
        .len()
        .filter(pl.col("len") == pl.col("len").max().over(["query_name", "strand"]))
        .select(["query_name", "strand", "#reference_name"])
    )

    df_minmax = (
        df.join(df_concensus_mapping, on=["query_name", "strand"], how="left")
        .rename(
            {
                "#reference_name_right": "ref_name",
            }
        )
        .group_by(["query_name", "ref_name", "strand"])
        .agg([pl.col("query_start").min(), pl.col("query_end").max()])
        .with_columns(query_len=pl.col("query_end") - pl.col("query_start"))
        .filter(pl.col("query_len") > args.len_thr)
        .sort("query_name")
        .select(
            [
                "query_name",
                "query_start",
                "query_end",
                "ref_name",
                "strand",
                "query_len",
            ]
        )
    )
    df_minmax.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    main()
