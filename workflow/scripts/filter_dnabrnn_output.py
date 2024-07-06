import re
import sys
import json
import argparse
import polars as pl
from typing import Any


# Reverse to not greedily match on chr2 instead of chr21
CHRS = [*[f"chr{i}" for i in range(22, 0, -1)], "chrX", "chrY"]
CHRS_ASAT_SEP = {"chr3", "chr4", "chr5", "chr7"}
CHRS_13_21 = {"chr13", "chr21"}
ORTS = ["fwd", "rev"]
RGX_CHR = re.compile("|".join(c for c in CHRS))
DEF_MERGE_REPEAT_DST_THR = 100_000
DEF_REPEAT_LEN_THR = [1_000, None]
DEF_REPEAT_TYPE = 2


def build_interval_expr(interval: tuple[int, int | None]) -> pl.Expr | None:
    try:
        thr_min, thr_max = interval
    except ValueError:
        return None

    if thr_min and thr_max:
        return (pl.col("rlen") > thr_min) & (pl.col("rlen") < thr_max)
    elif thr_min:
        return pl.col("rlen") > thr_min
    else:
        return pl.col("rlen") < thr_max


def main():
    ap = argparse.ArgumentParser(
        description="Filter alpha-satellite repeats from dna-brnn output."
    )
    ap.add_argument("-i", "--input", required=True, help="Input dna-brnn output.")
    ap.add_argument(
        "-o",
        "--output",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="Output filtered dna-brnn output.",
    )
    ap.add_argument(
        "-t", "--thresholds", help="Input JSON file with repeat length thresholds."
    )
    ap.add_argument(
        "-c",
        "--chr",
        choices=CHRS,
        required=True,
        help="Chromosome to filter for.",
        metavar="{chr1|...|chr22|chrX|chrY}",
    )
    ap.add_argument(
        "-r",
        "--repeat_type",
        type=int,
        default=DEF_REPEAT_TYPE,
        help="Repeat type to filter for.",
    )

    args = ap.parse_args()

    try:
        df = pl.read_csv(
            args.input,
            separator="\t",
            has_header=False,
            new_columns=["ctg", "start", "end", "rtype"],
        )
    except pl.exceptions.NoDataError:
        return

    if args.thresholds:
        thresholds = json.load(open(args.thresholds))
    else:
        thresholds = {}

    selected_chr = args.chr
    selected_repeat_type = args.repeat_type
    repeat_len_thresholds = thresholds.get("repeat_len_thr", {})
    default_repeat_len_threshold = tuple(
        thresholds.get("default_repeat_len_thr", DEF_REPEAT_LEN_THR)
    )
    merge_repeat_dst_thr = thresholds.get("merge_repeat_dst_thr", {})
    default_merge_repeat_dst_thr = thresholds.get(
        "default_merge_repeat_dst_thr", DEF_MERGE_REPEAT_DST_THR
    )

    # Parse contig start and stop coords.
    # Only one ':' should be in full contig name.
    # ex. sample_chr_ctgname:start-end
    lf = (
        df.lazy()
        .with_columns(
            ctg=pl.col("ctg").str.split_exact(by=":", n=1),
            rlen=pl.col("end") - pl.col("start"),
        )
        .unnest("ctg")
        .rename({"field_0": "ctg", "field_1": "coords"})
        .with_columns(pl.col("coords").str.split_exact(by="-", n=1))
        .unnest("coords")
        .rename({"field_0": "ctg_start", "field_1": "ctg_end"})
        .cast({"ctg_start": pl.Int64, "ctg_end": pl.Int64})
    )
    # Adjust start and end positions with contig full coords.
    lf = lf.with_columns(
        start=pl.col("start") + pl.col("ctg_start"),
        end=pl.col("end") + pl.col("ctg_start"),
    )

    df = lf.select("ctg", "start", "end", "rtype", "rlen").collect()
    dfs = []

    for df_ctg_name, df_ctg in df.group_by(["ctg"]):
        df_ctg_name = df_ctg_name[0]
        chr_len_thr_exprs = []

        if mtch_chr_name := re.search(RGX_CHR, df_ctg_name):
            chr_name = mtch_chr_name.group().strip("_")
            if chr_name != selected_chr:
                continue

            # Build repeat len filtering expressions.
            for len_thr in repeat_len_thresholds.get(
                chr_name, [default_repeat_len_threshold]
            ):
                thr_expr = build_interval_expr(tuple(len_thr))
                if isinstance(thr_expr, pl.Expr):
                    chr_len_thr_exprs.append(thr_expr)
        else:
            continue

        # Apply repeat len filters.
        df_ctg_repeats = (
            df_ctg.filter(pl.any_horizontal(chr_len_thr_exprs))
            .filter(pl.col("rtype") == selected_repeat_type)
            .select("ctg", "start", "end", "rtype", "rlen")
        )

        # Merge adjacent rows within some dst
        all_rows: list[dict[str, Any]] = list(df_ctg_repeats.iter_rows(named=True))
        curr_pos = 0
        chr_merge_repeat_dst_thr = merge_repeat_dst_thr.get(
            chr_name, default_merge_repeat_dst_thr
        )
        while True:
            try:
                curr_row = all_rows[curr_pos]
                next_row = all_rows[curr_pos + 1]
                dst = next_row["start"] - curr_row["end"]
                if dst < chr_merge_repeat_dst_thr:
                    new_row = {
                        "ctg": curr_row["ctg"],
                        "start": curr_row["start"],
                        "end": next_row["end"],
                        "rtype": curr_row["rtype"],
                        "rlen": next_row["end"] - curr_row["start"],
                    }
                    all_rows.pop(curr_pos + 1)
                    all_rows.pop(curr_pos)
                    all_rows.insert(curr_pos, new_row)
                else:
                    curr_pos += 1
            except IndexError:
                break

        df_ctg_compressed_repeats = pl.DataFrame(
            all_rows, schema=["ctg", "start", "end", "rtype", "rlen"]
        )
        if df_ctg_compressed_repeats.is_empty():
            continue

        # If a chromosome separated by hsat/other repeat or has multiple arrays, set repeat bounds differently.
        # * asat sep - all repeats filtering smaller alpha-sat repeats.
        # * chr13/21 - largest repeat and 2 adjacent rows.
        # * otherwise - largest repeat. single hor array.
        if chr_name in CHRS_ASAT_SEP:
            df_ctg_compressed_repeats = df_ctg_compressed_repeats.filter(
                pl.col("rlen") >= 50_000
            )
        elif chr_name in CHRS_13_21:
            # Assumes that already correctly oriented!
            # 2 repeats away for chr13/21
            largest_rlen_row_num = df_ctg_compressed_repeats.with_row_index().filter(
                pl.col("rlen") == pl.col("rlen").max()
            )["index"][0]
            df_ctg_compressed_repeats = (
                df_ctg_compressed_repeats.with_row_index()
                .filter(
                    (pl.col("index") >= largest_rlen_row_num - 2)
                    & (pl.col("index") <= largest_rlen_row_num)
                )
                .drop("index")
            )
        else:
            df_ctg_compressed_repeats = df_ctg_compressed_repeats.filter(
                pl.col("rlen") == pl.col("rlen").max()
            )
        dfs.append(df_ctg_compressed_repeats)

    if dfs:
        df = pl.concat(dfs)
    else:
        # Allow creating an empty dataframe.
        df = pl.DataFrame()

    df.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())
