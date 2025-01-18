import sys
import argparse
import polars as pl

OUTPUT_FIELDS = [
    "new_name",
    "new_cen_st",
    "new_cen_end",
    "is_partial",
    "ctg_name",
    "ctg_len",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i", "--statuses", nargs="+", help="All cen status files by sample."
    )
    ap.add_argument("-f", "--faidx", required=True, help="Sample fasta index.")
    ap.add_argument(
        "-o",
        "--output",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help=f"Output rename key. Format: {OUTPUT_FIELDS}",
    )
    args = ap.parse_args()

    # join_cen_status_and_og_status
    df_all_statuses = pl.concat(
        pl.read_csv(
            file,
            separator="\t",
            new_columns=["old_name", "new_name", "ort", "is_partial"],
            has_header=False,
        )
        for file in args.statuses
    )

    df_statuses = (
        df_all_statuses
        # No underscores or |. These are reserved characters for naming.
        # We expect these fields are splitting.
        .with_columns(
            name_elems=pl.col("old_name")
            .str.extract_groups(r"(.*?)_(chr[0-9XY]+)_(.*?):(.*?)$")
            .struct.rename_fields(["sample", "chrom", "ctg_name", "coords"]),
        )
        .unnest("name_elems")
        .with_columns(
            pl.col("coords")
            .str.split_exact(by="-", n=1)
            .struct.rename_fields(names=["cen_st", "cen_end"])
        )
        .unnest("coords")
        .cast({"cen_st": pl.Int64, "cen_end": pl.Int64})
    )

    df_fai = pl.read_csv(
        args.faidx,
        has_header=False,
        new_columns=["ctg_name", "ctg_len", "a", "b", "c"],
        separator="\t",
    ).drop("a", "b", "c")

    # get_final_rename_key
    df_status_w_ctg_len = df_statuses.join(df_fai, on=["ctg_name"]).with_columns(
        new_name=pl.when(pl.col("ort") == "rev")
        .then(pl.col("sample") + "_rc-" + pl.col("chrom") + "_" + pl.col("ctg_name"))
        .otherwise(pl.col("sample") + "_" + pl.col("chrom") + "_" + pl.col("ctg_name")),
        new_cen_st=pl.when(pl.col("ort") == "rev")
        .then(pl.col("ctg_len") - pl.col("cen_end") + 1)
        .otherwise(pl.col("cen_st")),
        new_cen_end=pl.when(pl.col("ort") == "rev")
        .then(pl.col("ctg_len") - pl.col("cen_st") + 1)
        .otherwise(pl.col("cen_end")),
    )
    df_status_w_ctg_len.select(OUTPUT_FIELDS).write_csv(
        args.output, include_header=False, separator="\t"
    )


if __name__ == "__main__":
    raise SystemExit(main())
