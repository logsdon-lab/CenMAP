import sys
import argparse
import polars as pl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default=sys.stdin, type=argparse.FileType("rb"))
    ap.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))
    ap.add_argument("--color_alr_alpha", type=str, default="#8B008B")
    ap.add_argument("--color_other", type=str, default="gray")
    ap.add_argument("--color_none", type=str, default="yellow")

    args = ap.parse_args()
    df = pl.read_csv(
        args.input,
        new_columns=["name", "start", "stop", "rtype", "action"],
        has_header=False,
        separator="\t",
    )

    dfs = []
    for name, df_grp in df.group_by(["name"]):
        name = name[0]
        df_grp = (
            df_grp.with_columns(
                pl.when(pl.col("rtype") == "ALR/Alpha")
                .then(pl.col("rtype"))
                .otherwise(pl.lit("Other"))
            )
            .with_columns(
                rle_id=pl.col("rtype").rle_id(),
                # Set color for other and plot.
                action=(
                    pl.when(pl.col("rtype") == "Other")
                    .then(pl.lit(f"plot:{args.color_other},ignore:absolute"))
                    .otherwise(pl.lit(f"plot:{args.color_alr_alpha}"))
                ),
            )
            .group_by(["rle_id"])
            .agg(
                pl.col("name").first(),
                pl.col("start").min(),
                pl.col("stop").max(),
                pl.col("rtype").first(),
                pl.col("action").first(),
            )
            .with_columns(len=pl.col("stop") - pl.col("start"))
            .sort(by="rle_id")
            # Handle edge cases where small ignored repeats between ALR. Allow misassemblies.
            .with_columns(
                pl.when(
                    (pl.col("rtype").shift(1) == "ALR/Alpha")
                    & (pl.col("rtype").shift(-1) == "ALR/Alpha")
                    & (pl.col("len") < 10_000)
                )
                .then(
                    pl.col("action").str.replace(",ignore:absolute", "", literal=True)
                )
                .otherwise(pl.col("action"))
            )
            .drop("rle_id", "len")
        )
        df_no_annotation = (
            df_grp.with_columns(diff=pl.col("start").shift(-1) - pl.col("stop"))
            .select(
                name=pl.lit(name),
                # Any set of rows where diff from first and second row is greater than 1 bp.
                start=pl.when(pl.col("diff") > 1)
                .then(pl.col("stop") + 1)
                .otherwise(None),
                stop=pl.when(pl.col("diff") > 1)
                .then(pl.col("start").shift(-1) - 1)
                .otherwise(None),
                rtype=pl.lit("None"),
                action=pl.lit("plot:yellow"),
            )
            .drop_nulls()
        )
        dfs.append(pl.concat([df_grp, df_no_annotation]).sort(by="start"))

    df_group_merged = pl.concat(dfs)
    df_group_merged.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())
