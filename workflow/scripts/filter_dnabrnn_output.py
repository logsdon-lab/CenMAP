import re
import sys
import json
import argparse
import polars as pl


CHRS = [*[f"chr{i}" for i in range(1, 23)], "chrX", "chrY"]
ORTS = ["fwd", "rev"]
RGX_CHR = re.compile("|".join(f"{c}_" for c in CHRS))
DEF_REPEAT_PERC_THR = 0.95
DEF_REPEAT_LEN_THR = [1_000, None]
DEF_REPEAT_TYPE = 2
DEF_DST_FROM_LARGEST = 4_000_000

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
        "-d",
        "--dst_from_largest",
        type=int,
        default=DEF_DST_FROM_LARGEST,
        help="Distance from largest repeat allowed.",
    )

    args = ap.parse_args()

    df = pl.read_csv(
        # "results/dna_brnn/K1463_2212/K1463_2212_centromeric_regions.renamed.fwd.bed",
        args.input,
        separator="\t",
        has_header=False,
        new_columns=["ctg", "start", "end", "rtype"],
    )
    if args.thresholds:
        # thresholds = json.load(open("config/dnabrnn_thresholds.json"))
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
    min_repeat_len_perc_threshold = thresholds.get(
        "min_repeat_len_perc", DEF_REPEAT_PERC_THR
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

        # Calculate percentile of repeats by repeat length.
        # Only allow repeats in the x-th percentile or above. This removes small repeats that would overextend the HOR array region. Only done if not all alpha-satellite.
        # Apply additional static repeat len filters.
        df_ctg = (
            df_ctg.filter(pl.any_horizontal(chr_len_thr_exprs))
            .with_columns(perc=pl.col("rlen").rank() / pl.col("rlen").len())
            .filter(
                pl.when((pl.col("rtype") == 2).all())
                .then(pl.lit(True))
                .otherwise(pl.col("perc") > min_repeat_len_perc_threshold)
                & (pl.col("rtype") == selected_repeat_type)
            )
            .select("ctg", "start", "end", "rtype", "rlen")
        )
        largest_row = df_ctg.filter(pl.col("rlen") == pl.col("rlen").max()).to_dict()

        # Allow only repeats some number of bps from the largest detected alpha-sat repeat.
        df_ctg = (
            df_ctg
            .filter(
                (pl.col("start") > largest_row["start"] - args.dst_from_largest) &
                (pl.col("end") < largest_row["end"] + args.dst_from_largest)
            )
        )

        dfs.append(df_ctg)

    if dfs:
        df = pl.concat(dfs)
    else:
        # Allow creating an empty dataframe.
        df = pl.DataFrame()

    df.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())
