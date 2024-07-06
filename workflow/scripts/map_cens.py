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
ACRO_CHRS = {"chr21", "chr22", "chr13", "chr14", "chr15"}


def main():
    ap = argparse.ArgumentParser("Map cens from alignments.")
    ap.add_argument(
        "-i",
        "--input",
        help="Input alignment bed intersected with chr p and q arms.",
        default=sys.stdin,
        type=argparse.FileType("rb"),
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
        pl.when(~pl.col("reference_name").is_in(ACRO_CHRS))
        .then(
            pl.all_horizontal(
                (pl.col("arm") == "p-arm").any(), (pl.col("arm") == "q-arm").any()
            ).over("query_name")
        )
        .otherwise(
            pl.all_horizontal((pl.col("arm") == "q-arm").any()).over("query_name")
        )
    )

    df_qarms = df.filter(pl.col("arm") == "q-arm").filter(
        pl.col("matches") == pl.col("matches").max().over(["query_name"])
    )
    df_concensus_mapping = (
        # Default to picking reference by number of matches
        df.filter(pl.col("matches") == pl.col("matches").max().over(["query_name"]))
        .join(df_qarms, on=["query_name"], how="left")
        .select("query_name", "reference_name", "reference_name_right", "arm_right")
        .unique()
        # But if has alignment to q-arm, take that instead.
        .with_columns(
            reference_name=pl.when(pl.col("arm_right") == "q-arm")
            .then(pl.col("reference_name_right"))
            .otherwise(pl.col("reference_name"))
        )
        .select(["query_name", "reference_name_right"])
    )
    df_ctg_groups = (
        df.join(df_concensus_mapping, on=["query_name"], how="left")
        .rename(
            {
                "reference_name_right": "final_reference_name",
            }
        )
        .group_by(["query_name", "final_reference_name"])
    )

    rows_ctg_grps: list[tuple[str, int, int, int, str]] = []
    for (qname, rname), df_ctg_grp in df_ctg_groups:
        ctg_start = df_ctg_grp["query_start"].min()
        ctg_end = df_ctg_grp["query_end"].max()
        ctg_len = ctg_end - ctg_start

        if ctg_len < args.len_thr:
            continue

        # Find arm mapping with highest number of matches.
        df_ctg_pqarm_mapping = df_ctg_grp.filter(pl.col("arm") != ".").filter(
            pl.col("matches") == pl.col("matches").max().over(["arm"])
        )

        # Reverse complement if needed.
        df_ort_check = df_ctg_pqarm_mapping.with_columns(
            exp_arm=pl.when(pl.col("query_start") == pl.col("query_start").min())
            .then(pl.lit("p-arm"))
            .otherwise(pl.lit("q-arm"))
        )
        if df_ort_check.is_empty():
            continue
        elif df_ort_check.shape[0] == 1:
            # Unable to determine without both arms. Assume already correctly oriented.
            is_not_rc = True
        else:
            is_not_rc = (df_ort_check["arm"] == df_ort_check["exp_arm"]).all()

        if is_not_rc:
            adj_ctg_start = ctg_start
            adj_ctg_end = ctg_end
        else:
            # Adjust coordinates if contig map in reverse ort.
            adj_ctg_start = df_ctg_grp["query_length"][0] - ctg_end
            adj_ctg_end = df_ctg_grp["query_length"][0] - ctg_start

        rows_ctg_grps.append(
            (
                qname,
                adj_ctg_start,
                adj_ctg_end,
                rname,
                adj_ctg_end - adj_ctg_start,
                not is_not_rc,
            )
        )

    df_minmax = pl.DataFrame(
        rows_ctg_grps,
        schema={
            "query_name": pl.String,
            "query_start": pl.Int64,
            "query_end": pl.Int64,
            "reference_name": pl.String,
            "query_len": pl.Int64,
            "reverse_complement": pl.Boolean,
        },
    )
    df_minmax.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    main()
