import re
import sys
import json
import argparse
import polars as pl
from typing import Any


CHRS = [*[f"chr{i}" for i in range(1, 23)], "chrX", "chrY"]
CHRS_MULTI_HOR_ARR = {"chr3", "chr4"}
ORTS = ["fwd", "rev"]
RGX_CHR = re.compile("|".join(f"{c}_" for c in CHRS))
DEF_MERGE_REPEAT_DST_THR = 10_000
DEF_REPEAT_LEN_THR = [1_000, None]
DEF_REPEAT_TYPE = 2
DEF_EDGES_BP = 500_000


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
        "--orientation",
        default="fwd",
        choices=ORTS,
        help="Orientation of contigs. Alters start and end calculations.",
    )
    ap.add_argument(
        "-c",
        "--chr",
        required=True,
        choices=CHRS,
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
    ap.add_argument(
        "-e",
        "--edges_bp",
        type=int,
        default=DEF_EDGES_BP,
        help="Base pairs required on edges. Filters out partial contigs.",
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

    ort = args.orientation
    selected_chr = args.chr
    selected_repeat_type = args.repeat_type
    repeat_len_thresholds = thresholds.get("repeat_len_thr", {})
    default_repeat_len_threshold = tuple(
        thresholds.get("default_repeat_len_thr", DEF_REPEAT_LEN_THR)
    )
    merge_repeat_dst_thr = thresholds.get(
        "merge_repeat_dst_thr", DEF_MERGE_REPEAT_DST_THR
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
    if ort == "rev":
        lf = lf.with_columns(
            start=pl.col("ctg_end") - pl.col("end"),
            end=pl.col("ctg_end") - pl.col("start"),
        )
    else:
        lf = lf.with_columns(
            start=pl.col("start") + pl.col("ctg_start"),
            end=pl.col("end") + pl.col("ctg_start"),
        )

    df = lf.select(
        "ctg", "start", "end", "rtype", "rlen", "ctg_start", "ctg_end"
    ).collect()
    dfs = []

    for df_ctg_name, df_ctg in df.group_by(["ctg"]):
        df_ctg_name = df_ctg_name[0]
        chr_len_thr_exprs = []

        ctg_start, ctg_end = df_ctg.row(0)[5], df_ctg.row(0)[6]

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
        rows: list[dict[str, Any]] = []
        for curr_row in df_ctg_repeats.iter_rows(named=True):
            try:
                prev_row = rows.pop()
                dst = curr_row["start"] - prev_row["end"]
                if dst < merge_repeat_dst_thr:
                    new_row = {
                        "ctg": curr_row["ctg"],
                        "start": prev_row["start"],
                        "end": curr_row["end"],
                        "rtype": curr_row["rtype"],
                        "rlen": curr_row["end"] - prev_row["start"],
                    }
                    rows.append(new_row)
                else:
                    rows.append(prev_row)
                    rows.append(curr_row)
            except IndexError:
                rows.append(curr_row)

        df_ctg_compressed_repeats = pl.DataFrame(
            rows, schema=["ctg", "start", "end", "rtype", "rlen"]
        )
        if df_ctg_compressed_repeats.is_empty():
            continue

        # If a chromosome that has multiple HOR arrays, set repeat bounds differently.
        # * multiple - use min and max coordinate of all repeats.
        # * single - use min and max coordinate of largest repeat.
        if chr_name in CHRS_MULTI_HOR_ARR:
            df_repeat_bounds = df_ctg_compressed_repeats.filter(
                pl.col("rlen") >= 100_000
            )
        else:
            df_repeat_bounds = df_ctg_compressed_repeats.filter(
                pl.col("rlen") == pl.col("rlen").max()
            )
        repeats_edge_start = df_repeat_bounds.row(0)[1] - DEF_EDGES_BP
        repeats_edge_end = df_repeat_bounds.row(-1)[2] + DEF_EDGES_BP

        # Ensure that repeats have at least some number of bps(default: 50 kbp) on either side.
        # TODO: Add test case to check if this works.
        if repeats_edge_start < ctg_start or repeats_edge_start < 0:
            continue
        elif repeats_edge_end > ctg_end:
            continue

        # Allow only repeats some number of bps from the largest detected alpha-sat repeat.
        try:
            df_ctg_compressed_repeats = df_ctg_compressed_repeats.filter(
                (pl.col("start") > repeats_edge_start)
                & (pl.col("end") < repeats_edge_end)
            )
        except IndexError:
            # If empty.
            pass
        dfs.append(df_ctg_compressed_repeats)

    if dfs:
        df = pl.concat(dfs)
    else:
        # Allow creating an empty dataframe.
        df = pl.DataFrame()

    df.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())
