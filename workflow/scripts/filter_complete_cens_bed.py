import sys
import argparse
import polars as pl


INFILE_COLS = ("name", "start", "end", "length")


def main():
    ap = argparse.ArgumentParser("Filter correct centromere bed from provided lists.")
    ap.add_argument(
        "-i",
        "--infile",
        type=str,
        help=f"Input centromere bed file with cols: {INFILE_COLS}",
    )
    ap.add_argument(
        "-l",
        "--lengths",
        type=str,
        help="Contig lengths. Contig names should match infile name col. Requires at least two columns. ex. Fasta index file",
    )
    ap.add_argument(
        "-c",
        "--correct_list",
        type=argparse.FileType("rt"),
        help="Correct cens list. Kept in bed.",
    )
    ap.add_argument(
        "-p",
        "--partial_list",
        type=argparse.FileType("rt"),
        help="Partial cens list. Removed from bed.",
    )
    ap.add_argument(
        "-r",
        "--reverse_key",
        type=argparse.FileType("rt"),
        help="Reverse cens key-value pairs. Kept in bed, replacement used, and coordinates adjusted.",
    )
    ap.add_argument(
        "-o",
        "--outfile",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="Output filtered and formatted bedfile.",
    )

    args = ap.parse_args()

    try:
        df_bed = pl.read_csv(
            args.infile, new_columns=INFILE_COLS, has_header=False, separator="\t"
        )
    except pl.exceptions.NoDataError:
        return

    df_fai = (
        pl.read_csv(args.lengths, has_header=False, separator="\t")
        .rename({"column_1": "name", "column_2": "length"})
        .select("name", "length")
    )
    partials = {line.strip() for line in args.partial_list.readlines()}
    correct = {line.strip() for line in args.correct_list.readlines()}
    reverse = dict(line.strip().split("\t") for line in args.reverse_key.readlines())

    df_final_bed = (
        df_bed.filter(
            (~pl.col("name").is_in(partials)) | (pl.col("name").is_in(correct))
        )
        .join(df_fai, how="left", on="name")
        .with_columns(
            new_name=pl.col("name").replace(reverse),
            new_start=pl.when(pl.col("name").is_in(reverse.keys()))
            .then(pl.col("length_right") - pl.col("end"))
            .otherwise(pl.col("start")),
            new_end=pl.when(pl.col("name").is_in(reverse.keys()))
            .then(pl.col("length_right") - pl.col("start"))
            .otherwise(pl.col("end")),
        )
        .select("new_name", "new_start", "new_end", "length")
    )
    df_final_bed.write_csv(args.outfile, include_header=False, separator="\t")


if __name__ == "__main__":
    raise SystemExit(main())
