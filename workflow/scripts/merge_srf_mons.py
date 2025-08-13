import sys
import argparse
import polars as pl


COLS = ("chrom", "st", "end", "name", "score", "strand", "tst", "tend", "item_rgb")


def main():
    ap = argparse.ArgumentParser(
        description="Merge srf stringdecomposed monomers into stv"
    )
    ap.add_argument(
        "-i",
        "--infile",
        type=argparse.FileType("rb"),
        required=True,
        help="Input stringdecomposer bed.",
    )
    ap.add_argument(
        "-o",
        "--outfile",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="Output stv bed.",
    )
    ap.add_argument("-t", "--thr_ident", help="Threshold identity.", default=70.0)

    args = ap.parse_args()
    (
        pl.read_csv(args.infile, separator="\t", has_header=False, new_columns=COLS)
        .filter(pl.col("score") > args.thr_ident)
        .with_columns(
            pl.col("chrom").str.extract(r"^(.+):"),
            chrom_st=pl.col("chrom")
            .str.extract(r":(\d+)-")
            .fill_null("0")
            .cast(pl.Int64),
        )
        .with_columns(
            pl.col("st") + pl.col("chrom_st"), pl.col("end") + pl.col("chrom_st")
        )
        .drop("chrom_st")
    )


if __name__ == "__main__":
    raise SystemExit(main())
