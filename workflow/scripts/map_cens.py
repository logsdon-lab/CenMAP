import sys
import argparse
import polars as pl

ALN_HEADER = [
    "reference_name",
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
    "mon_reference_name",
    "mon_start",
    "mon_end",
    "arm",
]


def main():
    ap = argparse.ArgumentParser("Map cens from alignments.")
    ap.add_argument(
        "-i",
        "--input",
        help="Input alignment bed intersected with chr p and q arms.",
        type=str,
    )
    ap.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))
    ap.add_argument("-t", "--len_thr", default=1_000_000, type=int)

    args = ap.parse_args()

    df = pl.read_csv(
        args.input, separator="\t", new_columns=ALN_HEADER, has_header=False
    )

    # To pass, a given query sequence must have both p and q arms mapped.
    # The acrocentrics are exceptions and only require the q-arm.
    df = df.filter(
        pl.when(
            ~pl.col("reference_name").is_in(
                ["chr21", "chr22", "chr13", "chr14", "chr15"]
            )
        )
        .then(
            pl.all_horizontal(
                (pl.col("arm") == "p-arm").any(), (pl.col("arm") == "q-arm").any()
            ).over("query_name")
        )
        .otherwise(
            pl.all_horizontal((pl.col("arm") == "q-arm").any()).over("query_name")
        )
    )

    # Choose concensus by taking highest number of matches .
    df_concensus_mapping = df.filter(
        pl.col("matches") == pl.col("matches").max().over(["query_name"])
    ).select(["query_name", "reference_name"])

    df_minmax = (
        df.join(df_concensus_mapping, on=["query_name"], how="left")
        .rename(
            {
                "reference_name_right": "final_reference_name",
            }
        )
        .group_by(["query_name", "final_reference_name", "strand"])
        .agg([pl.col("query_start").min(), pl.col("query_end").max()])
        .with_columns(query_len=pl.col("query_end") - pl.col("query_start"))
        .filter(pl.col("query_len") > args.len_thr)
        .sort("query_name")
        .select(
            [
                "query_name",
                "query_start",
                "query_end",
                "final_reference_name",
                "strand",
                "query_len",
            ]
        )
    )
    df_minmax.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    main()
