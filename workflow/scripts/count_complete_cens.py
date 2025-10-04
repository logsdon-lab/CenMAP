import sys
import argparse

import polars as pl


DEF_IN_COLS = ["ctg", "start", "end", "length"]
DEF_OUT_SCHEMA = {"sample": pl.String, "len": pl.UInt32, "perc": pl.Float64}
DEF_CHRS = list(reversed([f"chr{i}" for i in [*range(1, 23), "X", "Y"]]))
DEF_N_CHR = (len(DEF_CHRS) * 2) - 2


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
    ap.add_argument(
        "-c",
        "--chroms",
        default=DEF_CHRS,
        type=str,
        nargs="+",
        help="Chromosome names.",
    )
    ap.add_argument(
        "-n",
        "--n_chroms",
        default=DEF_N_CHR,
        help="Number of chromosomes in diploid organism.",
    )

    args = ap.parse_args()
    chroms = args.chroms
    n_chroms = args.n_chroms
    no_chrom = "all" in args.chroms

    # Include rc- in pattern
    if not no_chrom:
        chroms.extend([f"rc-{chrom}" for chrom in chroms])
        rgx_chrom = "|".join([*chroms, "-"])
        rgx_name_groups = (
            r"^(?<sample>.*?)_(?<chr>(" + rgx_chrom + r")*)_(?<contig_name>.*?)$"
        )
    else:
        rgx_name_groups = r"^(?<sample>.*?)_(?<contig_name>.*?)$"
    df = (
        pl.read_csv(
            args.input,
            separator="\t",
            has_header=False,
            new_columns=DEF_IN_COLS,
        )
        .group_by("ctg")
        .agg(pl.sum("length").alias("length"))
        .with_columns(mtch_ctg=pl.col("ctg").str.extract_groups(rgx_name_groups))
        .unnest("mtch_ctg")
    )

    df_cnts = (
        df.group_by("sample")
        .len()
        .with_columns(perc=((pl.col("len") / n_chroms) * 100).round(1))
    )
    df_cnts = pl.concat(
        [
            df_cnts,
            pl.DataFrame(
                data={
                    "sample": "all",
                    "len": df_cnts["len"].sum(),
                    "perc": round(
                        (df_cnts["len"].sum() / (df_cnts.shape[0] * n_chroms)) * 100, 1
                    ),
                },
                schema=DEF_OUT_SCHEMA,
            ),
        ]
    ).select(DEF_OUT_SCHEMA.keys())
    df_cnts.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())
