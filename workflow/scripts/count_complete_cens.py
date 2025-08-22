import sys
import argparse

import polars as pl


DEF_IN_COLS = ["ctg", "start", "end", "length"]
DEF_OUT_SCHEMA = {"sample": pl.String, "len": pl.UInt32, "perc": pl.Float64}
NUM_CHRS = 46


def main():
    ap = argparse.ArgumentParser(
        description="Count complete centromeres from HOR array length."
    )
    ap.add_argument(
        "-i",
        "--input",
        required=True,
        type=str,
        help=f"Input HOR array length by contig with sample, chromosome, and contig/haplotype in name (ex. HG00171_chr1_h2tg000057l#1-9773347:114672-7639070). Expects TSV with no header and the fields: {DEF_IN_COLS}",
    )
    ap.add_argument(
        "-o",
        "--output",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help=f"Output centromere counts as TSV with no header and the fields: {list(DEF_OUT_SCHEMA.keys())}",
    )

    args = ap.parse_args()
    df = (
        pl.read_csv(
            args.input,
            separator="\t",
            has_header=False,
            new_columns=DEF_IN_COLS,
        )
        .group_by("ctg")
        .agg(pl.sum("length").alias("length"))
        # TODO: Replace with argparse regex for chroms. Should be based on config chroms.
        .with_columns(
            pl.col("ctg").str.extract_groups(
                r"(.*?)_(rc-chr[XY0-9]+|chr[XY0-9]+)_(.*?)"
            )
        )
        .unnest("ctg")
        .rename({"1": "sample", "2": "chr", "3": "ctg"})
    )

    df_cnts = (
        df.group_by("sample")
        .len()
        .with_columns(perc=((pl.col("len") / NUM_CHRS) * 100).round(1))
    )
    df_cnts = pl.concat(
        [
            df_cnts,
            pl.DataFrame(
                data={
                    "sample": "all",
                    "len": df_cnts["len"].sum(),
                    "perc": round(
                        (df_cnts["len"].sum() / (df_cnts.shape[0] * NUM_CHRS)) * 100, 1
                    ),
                },
                schema=DEF_OUT_SCHEMA,
            ),
        ]
    )
    df_cnts.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())
