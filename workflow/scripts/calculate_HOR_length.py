import sys
import argparse
import polars as pl


def main():
    ap = argparse.ArgumentParser("Estimate HOR array length from HumAS-HMMER output.")
    ap.add_argument(
        "-i", "--input", help="Input bed file made from HumAS-HMMER output.", type=str
    )
    ap.add_argument(
        "-o",
        "--output",
        help="Output bed file with chr, start, stop, and len columns.",
        default=sys.stdout,
        type=argparse.FileType("wt"),
    )
    ap.add_argument(
        "--bp_jump_thr",
        help="Base pair jump threshold to group by",
        type=int,
        default=20000,
    )

    args = ap.parse_args()

    df = pl.read_csv(
        args.input,
        separator="\t",
        columns=[0, 1, 2, 3, 5],
        new_columns=["chr", "start", "stop", "hor", "strand"],
        has_header=False,
    )

    dfs = []
    for chr_name, df_chr in df.group_by(["chr"]):
        chr_name = chr_name[0]
        df_live_hor = df_chr.filter(pl.col("hor").str.contains("L"))

        df_bp_jumps = df_live_hor.with_columns(
            diff=pl.col("start") - pl.col("start").shift(1)
        ).filter(pl.col("diff") > args.bp_jump_thr)

        if df_bp_jumps.is_empty():
            dfs.append(
                df_live_hor.with_columns(
                    chr_name=pl.lit(chr_name),
                    start_pos=pl.col("start").min(),
                    stop_pos=pl.col("stop").min(),
                    len=pl.col("stop").max() - pl.col("start").min(),
                )
            )
            continue

        range_coords = []
        for i, row in enumerate(df_bp_jumps.iter_rows()):
            prev_row = pl.DataFrame() if i == 0 else df_bp_jumps.slice(i - 1)
            next_row = df_bp_jumps.slice(i + 1)

            if prev_row.is_empty():
                range_coords.append((df_chr.get_column("start").min(), row[1]))

            range_coords.append((row[1], row[2]))

            if next_row.is_empty():
                range_coords.append((row[2], df_chr.get_column("stop").max()))
            else:
                range_coords.append((row[2], next_row.get_column("start")[0]))

        breakpoint()


if __name__ == "__main__":
    raise SystemExit(main())
