import sys
import json
import argparse
import polars as pl

from enum import Enum, auto
from typing import Any


# Reverse to not greedily match on chr2 instead of chr21
CHRS = [*[f"chr{i}" for i in range(22, 0, -1)], "chrX", "chrY"]
CHRS_ASAT_SEP = {"chr3", "chr4", "chr5"}
CHRS_13_21 = {"chr13", "chr21"}
ORTS = ["fwd", "rev"]
DEF_ARR_LEN_THR = 30_000
DEF_MERGE_REPEAT_DST_THR = 100_000
DEF_REPEAT_LEN_THR = [1_000, None]
DEF_REPEAT_TYPE = 2
DEF_COMP_REPEAT_LEN_THR = 100_000


class MergeMode(Enum):
    """
    Limit bp that can be merged.
    """

    NoLimit = auto()
    Limit = auto()


class MergeOrt(Enum):
    Fwd = auto()
    Rev = auto()


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
    ap.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input dna-brnn output.",
        type=argparse.FileType("rb"),
    )
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
        help="Chromosome of sequence for dna-brnn output.",
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
    arr_len_thresholds = thresholds.get("arr_len_thr", {})
    def_arr_len_threshold = arr_len_thresholds.get("default", DEF_ARR_LEN_THR)

    repeat_len_thresholds = thresholds.get("repeat_len_thr", {})
    default_repeat_len_threshold = tuple(
        thresholds.get("default_repeat_len_thr", DEF_REPEAT_LEN_THR)
    )
    merge_repeat_dst_thr = thresholds.get("merge_repeat_dst_thr", {})
    default_merge_repeat_dst_thr = thresholds.get(
        "default_merge_repeat_dst_thr", DEF_MERGE_REPEAT_DST_THR
    )
    merge_mode = thresholds.get("merge_mode", {})
    default_merge_mode = thresholds.get("default_merge_mode", "NoLimit")
    merge_ort = thresholds.get("merge_ort", {})
    default_merge_ort = thresholds.get("default_merge_ort", "Fwd")
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

    for ctg_name, df_ctg in df.group_by(["ctg"]):
        ctg_name = ctg_name[0]
        chr_len_thr_exprs = []

        # Build repeat len filtering expressions.
        for len_thr in repeat_len_thresholds.get(
            selected_chr, [default_repeat_len_threshold]
        ):
            thr_expr = build_interval_expr(tuple(len_thr))
            if isinstance(thr_expr, pl.Expr):
                chr_len_thr_exprs.append(thr_expr)

        # Apply repeat len filters.
        df_ctg_repeats = (
            df_ctg.filter(pl.any_horizontal(chr_len_thr_exprs))
            .filter(pl.col("rtype") == selected_repeat_type)
            .select("ctg", "start", "end", "rtype", "rlen")
        )

        # Merge adjacent rows within some dst
        chr_merge_ort = merge_ort.get(selected_chr, default_merge_ort)
        if MergeOrt[chr_merge_ort] == MergeOrt.Rev:
            iter_rows = df_ctg_repeats.reverse().iter_rows(named=True)
        else:
            iter_rows = df_ctg_repeats.iter_rows(named=True)
        all_rows: list[dict[str, Any]] = list(iter_rows)
        curr_pos = 0
        chr_merge_repeat_dst_thr = merge_repeat_dst_thr.get(
            selected_chr, default_merge_repeat_dst_thr
        )
        chr_merge_mode = merge_mode.get(selected_chr, default_merge_mode)
        while True:
            try:
                curr_row = all_rows[curr_pos]
                next_row = all_rows[curr_pos + 1]
                if MergeOrt[chr_merge_ort] == MergeOrt.Fwd:
                    dst = next_row["start"] - curr_row["end"]
                    start, end = curr_row["start"], next_row["end"]
                    rlen = next_row["end"] - curr_row["start"]
                else:
                    dst = curr_row["start"] - next_row["end"]
                    start, end = next_row["start"], curr_row["end"]
                    rlen = curr_row["end"] - next_row["start"]

                if dst < chr_merge_repeat_dst_thr:
                    new_row = {
                        "ctg": curr_row["ctg"],
                        "start": start,
                        "end": end,
                        "rtype": curr_row["rtype"],
                        "rlen": rlen,
                    }
                    all_rows.pop(curr_pos + 1)
                    all_rows.pop(curr_pos)
                    all_rows.insert(curr_pos, new_row)

                    # Reduce distance that can be merged.
                    if MergeMode[chr_merge_mode] == MergeMode.Limit:
                        chr_merge_repeat_dst_thr -= dst
                else:
                    curr_pos += 1
            except IndexError:
                break

        df_ctg_compressed_repeats = pl.DataFrame(
            all_rows, schema=["ctg", "start", "end", "rtype", "rlen"]
        )
        if MergeOrt[chr_merge_ort] == MergeOrt.Rev:
            df_ctg_compressed_repeats = df_ctg_compressed_repeats.reverse()

        if df_ctg_compressed_repeats.is_empty():
            continue

        # If a chromosome separated by hsat/other repeat or has multiple arrays, set repeat bounds differently.
        # * asat sep - all repeats filtering smaller alpha-sat repeats.
        # * chr13/21 - largest repeat and 1 adjacent repeat.
        # * otherwise - largest repeat. single hor array.
        if selected_chr in CHRS_ASAT_SEP:
            df_ctg_compressed_repeats = df_ctg_compressed_repeats.filter(
                pl.col("rlen") >= DEF_COMP_REPEAT_LEN_THR
            )
        elif selected_chr in CHRS_13_21:
            df_ctg_compressed_repeats = df_ctg_compressed_repeats.filter(
                pl.col("rlen") >= DEF_COMP_REPEAT_LEN_THR
            ).with_row_index()
            df_largest_rlen_row = df_ctg_compressed_repeats.filter(
                pl.col("rlen") == pl.col("rlen").max()
            )
            try:
                largest_rlen_row_num = df_largest_rlen_row["index"][0]
            except IndexError:
                sys.stderr.write(
                    f"{ctg_name} doesn't have a alpha-satellite repeat larger than {DEF_COMP_REPEAT_LEN_THR:,} bp.\n"
                )
                continue

            # 1 repeats away for main HOR array chr13/21
            # Get side that trims the most.
            if df_ctg_compressed_repeats["index"].median() > largest_rlen_row_num:
                adj_row_num = largest_rlen_row_num + 1
            elif df_ctg_compressed_repeats["index"].median() < largest_rlen_row_num:
                adj_row_num = largest_rlen_row_num - 1
            else:
                adj_row_num = None

            if isinstance(adj_row_num, int):
                df_ctg_compressed_repeats = pl.concat(
                    [
                        df_largest_rlen_row,
                        # Get adjacent row and trim 500kbp ahead.
                        # This replaces closest row with boundaries which will be min/maxed in following scripts.
                        df_ctg_compressed_repeats.filter(
                            pl.col("index") == adj_row_num
                        ).with_columns(
                            start=pl.col("end") + 500_000,
                            end=pl.col("end") + 500_000,
                        ),
                    ]
                ).sort(by="start")
            df_ctg_compressed_repeats = df_ctg_compressed_repeats.drop("index")
        else:
            df_ctg_compressed_repeats = df_ctg_compressed_repeats.filter(
                pl.col("rlen") == pl.col("rlen").max()
            )
        # Aggregate and bedminmax all.
        # Select cols and calculate length.
        # Add 500 kbp on both ends.
        # Take only repeats greater than some value.
        # Take abs value.
        df_ctg_compressed_repeats = (
            df_ctg_compressed_repeats.group_by(["ctg", "rtype"])
            .agg(pl.col("start").min(), pl.col("end").max())
            .sort(by="start")
            .with_row_index()
            .with_columns(
                start=pl.when(pl.col("index") == 0)
                .then(pl.col("start") - 500_000)
                .otherwise(pl.col("start"))
                .clip(0, pl.col("start").max()),
                end=pl.when(pl.col("index") == df_ctg_compressed_repeats.shape[0])
                .then(pl.col("end") + 500_000)
                .otherwise(pl.col("end")),
                rlen=pl.col("end") - pl.col("start"),
            )
            .filter(
                pl.col("rlen")
                > arr_len_thresholds.get(selected_chr, def_arr_len_threshold)
            )
            .select("ctg", "start", "end", "rtype", "rlen")
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
