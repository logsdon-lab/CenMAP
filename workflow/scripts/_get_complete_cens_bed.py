import argparse
import polars as pl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--statuses", nargs="+", help="All cen status files by sample.")
    ap.add_argument("--rename_keys", nargs="+", help="All rename keys by sample.")
    ap.add_argument("--faidx", help="Sample fasta index.")
    ap.add_argument("--output", nargs="+")

    args = ap.parse_args()

    # join_cen_status_and_og_status
    df_statuses = pl.concat(
        pl.read_csv(
            file,
            separator="\t",
            new_columns=["old_name", "new_name", "ort", "is_partial"],
            has_header=False,
        )
        .with_columns(
            new_name=pl.col("old_name").str.extract(r"^(.*?):"),
            coords=pl.col("old_name").str.extract(r":(.*?)$"),
        )
        .select("new_name", "ort", "is_partial", "coords")
        for file in args.statuses
    )
    df_rename_key = pl.concat(
        pl.read_csv(
            file,
            separator="\t",
            new_columns=[
                "old_name",
                "new_name",
                "coords",
                "sample",
                "chrom",
                "is_reversed",
            ],
            has_header=False,
        )
        .with_columns(
            ort=pl.when(pl.col("is_reversed") == "false")
            .then(pl.lit("fwd"))
            .otherwise("rev")
        )
        .select("new_name", "ort")
        for file in args.rename_keys
    )

    df_status_rename_key = df_statuses.join(df_rename_key, on=["new_name"])
    print(df_status_rename_key)
    # get_final_rename_key

    # make_complete_cens_bed


if __name__ == "__main__":
    raise SystemExit(main())
