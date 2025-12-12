import sys
import argparse
import polars as pl
import scipy.signal
import intervaltree as it
from collections import defaultdict


def group_by_dst(df: pl.DataFrame, dst: int, group_name: str) -> pl.DataFrame:
    try:
        df = df.drop("index")
    except pl.exceptions.ColumnNotFoundError:
        pass
    return (
        df.with_columns(
            # c1  st1 (end1)
            # c1 (st2) end2
            dst_behind=(pl.col("st") - pl.col("end").shift(1)).fill_null(0),
            dst_ahead=(pl.col("st").shift(-1) - pl.col("end")).fill_null(0),
        )
        .with_row_index()
        .with_columns(
            **{
                # Group HOR units based on distance.
                group_name: pl.when(pl.col("dst_behind").le(dst))
                # We assign 0 if within merge dst.
                .then(pl.lit(0))
                # Otherwise, give unique index.
                .otherwise(pl.col("index") + 1)
                # Then create run-length ID to group on.
                # Contiguous rows within distance will be grouped together.
                .rle_id()
            },
        )
        .with_columns(
            # Adjust groups in scenarios where should join group ahead or behind but given unique group.
            # B:64617 A:52416 G:1
            # B:52416 A:1357  G:2 <- This should be group 3.
            # B:1357  A:1358  G:3
            pl.when(pl.col("dst_behind").le(dst) & pl.col("dst_ahead").le(dst))
            .then(pl.col(group_name))
            .when(pl.col("dst_behind").le(dst))
            .then(pl.col(group_name).shift(1))
            .when(pl.col("dst_ahead").le(dst))
            .then(pl.col(group_name).shift(-1))
            .otherwise(pl.col(group_name))
        )
    )


def main():
    ap = argparse.ArgumentParser(
        description="Determine if centromeres are valid based on kmer shannon index and ALR overlap. A complete centromere should have dip in entropy overlapping putative ALR."
    )
    ap.add_argument(
        "-i",
        "--infile",
        required=True,
        type=argparse.FileType("rb"),
        help="BED9 file with kmer shannon index as score.",
    )
    ap.add_argument(
        "-b",
        "--bed",
        type=argparse.FileType("rb"),
        help="BED4 annotation file with no header. Chromosome name does not contain contig coordinates. Absolute coordinates",
    )
    ap.add_argument(
        "-o",
        "--outfile",
        type=argparse.FileType("wt"),
        default=sys.stdout,
        help="BED3 of valid centromeres as coords,start,end. Extracted from name column of entropy bed.",
    )
    ap.add_argument(
        "-p",
        "--prop_valid",
        type=float,
        default=0.9,
        help="Proportion of dips in entropy over all putative ALR containing regions required to be valid.",
    )
    ap.add_argument(
        "-d",
        "--dst",
        type=int,
        default=1_000_000,
        help="Merge distance.",
    )
    args = ap.parse_args()

    try:
        df_entropy = (
            pl.read_csv(
                args.infile,
                separator="\t",
                has_header=False,
                new_columns=[
                    "chrom",
                    "st",
                    "end",
                    "name",
                    "score",
                    "strand",
                    "tst",
                    "tend",
                    "item_rgb",
                ],
            )
            .with_columns(
                mtches=pl.col("chrom").str.extract_groups(
                    r"^(?<chrom_name>.+):(?<chrom_st>\d+)-(?<chrom_end>\d+)$"
                )
            )
            .unnest("mtches")
            .cast({"chrom_st": pl.Int64, "chrom_end": pl.Int64})
            .with_columns(chrom_st=pl.col("chrom_st").fill_null(0))
            # Convert to absolute coordinates
            .with_columns(
                pl.col("st") + pl.col("chrom_st"),
                pl.col("end") + pl.col("chrom_st"),
                pl.col("tst") + pl.col("chrom_st"),
                pl.col("tend") + pl.col("chrom_st"),
                chrom=pl.col("chrom_name"),
            )
        )
    except pl.exceptions.NoDataError:
        return

    df_rm = pl.read_csv(
        args.bed,
        separator="\t",
        has_header=False,
        columns=[0, 1, 2, 3],
        new_columns=["chrom", "st", "end", "rtype"],
    )

    itree_rm = defaultdict(it.IntervalTree)
    for chrom, st, end, _ in df_rm.filter(pl.col("rtype") == "ALR/Alpha").iter_rows():
        itree_rm[chrom].add(it.Interval(st, end))

    df_entropy = (
        df_entropy.sort(by=["chrom", "st"])
        # Regions with any drop in entropy are candidates.
        # We'll intersect these with putative ALR regions (171 * n)
        .with_columns(
            score=pl.when(pl.col("score") < 1.0)
            .then(pl.lit(0.0))
            .otherwise(pl.lit(1.0))
        )
        # Get minimum and maximum position of zero regions.
        .with_columns(rle_score=pl.col("score").rle_id().over("chrom"))
        .group_by(["chrom", "rle_score", "chrom_st", "chrom_end"])
        .agg(
            pl.col("st").min(),
            pl.col("end").max(),
            pl.col("name").first(),
            pl.col("score").first(),
            pl.col("strand").first(),
            pl.col("tst").min(),
            pl.col("tend").max(),
            pl.col("item_rgb").first(),
        )
        .sort(by=["chrom", "st"])
        .drop("rle_score")
    )

    for chrom, df_chrom_entropy in df_entropy.group_by(
        ["chrom", "chrom_st", "chrom_end"]
    ):
        df_chrom_entropy = df_chrom_entropy.with_row_index()
        chrom, chrom_st, chrom_end = chrom
        # Find dips.
        dips, _ = scipy.signal.find_peaks(-df_chrom_entropy["score"].to_numpy())
        df_dips = df_chrom_entropy.filter(pl.col("index").is_in(dips))

        # Get total length overlapping ALR and has entropy not equal to 1.0 in window.
        total_length = 0
        for st, end, score in df_chrom_entropy.select("st", "end", "score").iter_rows():
            length = end - st
            # 1.0 is highest entropy. Anything below is allowed.
            if itree_rm[chrom].overlaps(st, end) and score != 1.0:
                total_length += length

        # Is dip and contains ALR
        valid_length = 0
        for st, end, score in df_dips.select("st", "end", "score").iter_rows():
            if itree_rm[chrom].overlaps(st, end) and score != 1.0:
                valid_length += end - st

        if total_length > 0 and ((valid_length / total_length) > args.prop_valid):
            print(f"{chrom}\t{chrom_st}\t{chrom_end}\t{chrom}:{chrom_st}-{chrom_end}")


if __name__ == "__main__":
    raise SystemExit(main())
