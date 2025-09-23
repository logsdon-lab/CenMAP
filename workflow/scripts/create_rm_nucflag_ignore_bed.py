import sys
import argparse
import polars as pl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, type=argparse.FileType("rb"))
    ap.add_argument("-r", "--regions", required=True, type=argparse.FileType("rb"))
    ap.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))

    args = ap.parse_args()
    df = pl.read_csv(
        args.input,
        new_columns=["name", "start", "stop", "rtype", "action"],
        has_header=False,
        separator="\t",
    ).sort(by=["name", "start"])
    df_regions = pl.read_csv(
        args.regions,
        new_columns=["name", "start", "stop"],
        columns=[0, 1, 2],
        has_header=False,
        separator="\t",
    )

    dfs = []
    for name, df_grp in df.group_by(["name"]):
        name = name[0]

        # Get whole region
        region_name, region_st, region_stop = df_regions.filter(
            pl.col("name") == name
        ).row(0)
        # Get annotated region
        annotated_st, annotated_stop = df_grp["start"].min(), df_grp["stop"].max()

        edge_rows = []
        if annotated_st > region_st:
            edge_rows.append(
                (region_name, region_st, annotated_st, "Other", "ignore:absolute")
            )
        if region_stop > annotated_stop:
            edge_rows.append(
                (region_name, annotated_stop, region_stop, "Other", "ignore:absolute")
            )

        # Ignore unannotated regions. This should be okay because scaffolded contigs will be removed by sh entropy step.
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
                    .then(pl.lit("ignore:absolute"))
                    .otherwise(None)
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
                .then(pl.col("action").str.replace("ignore:absolute", "", literal=True))
                .otherwise(pl.col("action"))
            )
            .filter(pl.col("action") != "")
            .drop("rle_id", "len")
        )
        df_unannotated = pl.DataFrame(edge_rows, schema=df_grp.schema, orient="row")
        dfs.append(pl.concat([df_grp, df_unannotated]).sort(by="start"))

    df_group_merged = pl.concat(dfs)
    df_group_merged.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())
