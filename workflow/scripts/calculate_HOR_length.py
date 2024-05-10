import re
import sys
import argparse
import polars as pl

RGX_CHR = re.compile(r"chr[0-9XY]{1,2}")


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
        default=100_000,
    )
    ap.add_argument(
        "--arr_len_thr",
        help="Length threshold to filter out",
        type=int,
        default=30_000,
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
    for ctg_name, df_chr in df.group_by(["chr"]):
        df_chr = df_chr.with_columns(len=pl.col("stop") - pl.col("start")).with_columns(
            mer=(pl.col("len") / 170).round()
        )
        ctg_name = ctg_name[0]
        chr_name = re.search(RGX_CHR, ctg_name).group()
        df_live_hor = df_chr.filter(pl.col("hor").str.contains("L"))

        # Specific edge case for chr8.
        if chr_name == "chr8" or chr_name == "chr10" or chr_name == "chr16":
            bp_jump_thr = 10_000
        elif chr_name == "chrY":
            bp_jump_thr = 2_000
        else:
            bp_jump_thr = args.bp_jump_thr

        df_bp_jumps = df_live_hor.with_columns(
            diff=pl.col("start") - pl.col("stop").shift(1)
        ).filter(pl.col("diff") > bp_jump_thr)

        if df_bp_jumps.is_empty():
            adj_start = df_live_hor.get_column("start").min()
            adj_stop = df_live_hor.get_column("stop").max()
            adj_len = adj_stop - adj_start

            if adj_len < args.arr_len_thr:
                continue

            dfs.append(
                pl.DataFrame(
                    {
                        "chr_name": ctg_name,
                        "start_pos": adj_start,
                        "stop_pos": adj_stop,
                        "len": adj_len,
                    }
                )
            )
            continue

        starts, stops = [], []
        for i, row in enumerate(df_bp_jumps.iter_rows()):
            prev_row = pl.DataFrame() if i == 0 else df_bp_jumps.slice(i - 1)
            next_row = df_bp_jumps.slice(i + 1)

            if prev_row.is_empty():
                starts.append(df_chr.get_column("start").min())
                stops.append(
                    df_chr.filter(pl.col("start") < row[1]).row(-1, named=True)["stop"]
                )

            if next_row.is_empty():
                starts.append(row[1])
                stops.append(df_chr.get_column("stop").max())
            else:
                starts.append(row[1])
                stops.append(
                    df_chr.filter(
                        pl.col("start") < next_row.get_column("start")[0]
                    ).row(-1, named=True)["stop"]
                )

        lens = []
        chr_mer_filter = None
        if chr_name == "chr10" or chr_name == "chr20":
            chr_mer_filter = pl.col("mer") >= 5
        elif chr_name == "chrY":
            chr_mer_filter = pl.col("mer") >= 30
        elif chr_name == "chr17":
            chr_mer_filter = pl.col("mer") >= 4

        for start, stop in zip(starts, stops):
            df_slice = (
                df_chr.filter(pl.col("start") >= start, pl.col("stop") <= stop)
                .with_columns(bp_jump=pl.col("start") - pl.col("stop").shift(1))
                .fill_null(0)
            )
            # Filter out mers based on chr.
            if chr_mer_filter is not None:
                df_slice = df_slice.filter(chr_mer_filter)

            if df_slice.is_empty():
                lens.append(0)
                continue
            df_slice_dst = (
                # df_slice.with_columns(len=pl.col("stop") - pl.col("start")).get_column("len").sum()
                df_slice.get_column("stop").max() - df_slice.get_column("start").min()
            )
            lens.append(df_slice_dst)

        lf = pl.LazyFrame(
            {
                "chr_name": ctg_name,
                "start_pos": starts,
                "stop_pos": stops,
                "len": lens,
            }
        )
        if (
            chr_name == "chr8"
            or chr_name == "chr10"
            or chr_name == "chr17"
            or chr_name == "chrY"
        ):
            arr_len_thr = 100_000
        else:
            arr_len_thr = args.arr_len_thr

        dfs.append(lf.filter(pl.col("len") > arr_len_thr).collect())

    df_all_dsts: pl.DataFrame = pl.concat(dfs)
    (
        df_all_dsts.with_columns(
            sort_idx=pl.col("chr_name")
            .str.extract("chr([0-9XY]+)")
            .replace({"X": "23", "Y": "24"})
            .cast(pl.Int32)
        )
        .sort(by="sort_idx")
        .select(["chr_name", "start_pos", "stop_pos", "len"])
        .write_csv(args.output, include_header=False, separator="\t")
    )


if __name__ == "__main__":
    raise SystemExit(main())
