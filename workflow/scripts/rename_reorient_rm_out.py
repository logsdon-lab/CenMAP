import sys
import argparse
import polars as pl

RM_FIELDS = (
    "idx",
    "div",
    "deldiv",
    "insdiv",
    "ctg",
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
)
KEY_FIELDS = (
    "new_name",
    "new_cen_st",
    "new_cen_end",
    "is_partial",
    "ctg_name",
    "ctg_len",
)


def main():
    ap = argparse.ArgumentParser(
        description="Adjust repeatmasker output by reorienting and adjust coordinates."
    )
    ap.add_argument(
        "-i", "--rm", required=True, help="Repeatmasker output with no header."
    )
    ap.add_argument(
        "-k",
        "--rename_key",
        help=f"Key used to rename. Format: {KEY_FIELDS}",
        required=True,
    )
    ap.add_argument(
        "-o",
        "--output",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="Output Repeatmasker output.",
    )
    args = ap.parse_args()

    df_rm = (
        pl.read_csv(
            args.rm,
            has_header=False,
            new_columns=RM_FIELDS,
            separator="\t",
            truncate_ragged_lines=True,
            schema_overrides={"right": pl.String},
        )
        .with_columns(
            pl.col("ctg")
            .str.split_exact(by=":", n=1)
            .struct.rename_fields(("ctg", "coords"))
        )
        .unnest("ctg")
        .with_columns(
            pl.col("coords")
            .str.split_exact(by="-", n=1)
            .struct.rename_fields(("cen_start", "cen_end"))
        )
        .unnest("coords")
        .cast({"cen_start": pl.Int64, "cen_end": pl.Int64})
    )
    df_rename_key = (
        pl.read_csv(
            args.rename_key, has_header=False, new_columns=KEY_FIELDS, separator="\t"
        )
        .with_columns(is_rc=(pl.col("new_name").str.contains("_rc-")))
        .with_columns(ctg=pl.col("new_name").str.replace("_rc-", "_"))
        .select("ctg", "new_name", "is_rc", "ctg_len")
    )

    df_final_rm = (
        df_rm.join(df_rename_key, on=["ctg"], how="inner")
        .with_columns(
            ctg=pl.when(pl.col("is_rc"))
            .then(pl.col("new_name"))
            .otherwise(pl.col("ctg")),
            cen_len=pl.col("cen_end") - pl.col("cen_start"),
        )
        .with_columns(
            start=pl.when(pl.col("is_rc"))
            .then(pl.col("cen_len") - pl.col("end") + 1)
            .otherwise(pl.col("start")),
            end=pl.when(pl.col("is_rc"))
            .then(pl.col("cen_len") - pl.col("start") + 1)
            .otherwise(pl.col("end")),
            ctg=pl.when(pl.col("is_rc"))
            .then(
                pl.col("ctg")
                + ":"
                + (pl.col("ctg_len") - pl.col("cen_end") + 1).cast(pl.String)
                + "-"
                + (pl.col("ctg_len") - pl.col("cen_start") + 1).cast(pl.String)
            )
            .otherwise(
                pl.col("ctg")
                + ":"
                + pl.col("cen_start").cast(pl.String)
                + "-"
                + pl.col("cen_end").cast(pl.String)
            ),
        )
        .drop("new_name", "is_rc", "cen_start", "cen_end", "cen_len", "ctg_len")
        .sort(by=["ctg", "start"])
    )
    df_final_rm.write_csv(args.output, include_header=False, separator="\t")


if __name__ == "__main__":
    raise SystemExit(main())
