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
        description="Determine if centromeres are valid based on shannon index and ALR overlap. A complete centromere should have dip in entropy overlapping ALR/Alpha."
    )
    ap.add_argument(
        "-i",
        "--infile",
        required=True,
        type=argparse.FileType("rb"),
        help="BED9 file with shannon index as score.",
    )
    ap.add_argument(
        "-r",
        "--repeatmasker",
        required=True,
        type=argparse.FileType("rb"),
        help="RepeatMasker output file with no header.",
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
        default=0.99,
        help="Proportion of dips in entropy over all ALR containing regions required to be valid. Majority by default.",
    )
    ap.add_argument(
        "-d",
        "--dst",
        type=int,
        default=1_000_000,
        help="Merge distance.",
    )
    ap.add_argument(
        "--trim_to_repeats",
        nargs="*",
        default=None,
        help="Trim coordinates to boundaries contain the largest block of these repeats.",
    )
    args = ap.parse_args()

    rm = args.repeatmasker
    # rm = "/project/logsdon_shared/projects/GIAB_HG008_TN/results/6-repeatmasker/repeats/HG008-N/HG008-N_rc-chr14_haplotype2-0000120:4461582-6946730.fa.out"
    bed = args.infile
    # bed = "/project/logsdon_shared/projects/GIAB_HG008_TN/exp_plot/HG008-N_rc-chr14_haplotype2-0000120:4461582-6946730.bed"
    thr_valid_prop = args.prop_valid

    try:
        df_entropy = pl.read_csv(
            bed,
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
    except pl.exceptions.NoDataError:
        return

    df_rm = pl.read_csv(
        rm,
        separator="\t",
        has_header=False,
        columns=[4, 5, 6, 9],
        new_columns=["chrom", "st", "end", "rtype"],
    )

    itree_rm = defaultdict(it.IntervalTree)
    for chrom, st, end, _ in df_rm.filter(pl.col("rtype") == "ALR/Alpha").iter_rows():
        itree_rm[chrom].add(it.Interval(st, end))

    df_entropy = (
        df_entropy.sort(by=["chrom", "st"])
        # Minmax regions of same entropy so peak calling easier.
        .with_columns(rle_score=pl.col("score").rle_id().over("chrom"))
        .group_by(["chrom", "rle_score"])
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
    )
    for chrom, df_chrom_entropy in df_entropy.group_by(["chrom"]):
        df_chrom_entropy = df_chrom_entropy.with_row_index()
        chrom: str = chrom[0]
        # Find dips.
        dips, _ = scipy.signal.find_peaks(-df_chrom_entropy["score"].to_numpy())
        df_dips = df_chrom_entropy.filter(pl.col("index").is_in(dips))

        # Get total length overlapping ALR and has zero entropy in window.
        total_length = 0
        for st, end, score in df_chrom_entropy.select("st", "end", "score").iter_rows():
            length = end - st
            if itree_rm[chrom].overlaps(st, end) and score == 0.0:
                total_length += length

        # Is dip and contains ALR
        # Handles case where a single dip due to LINE or ALU insertion. Will create small spike in entropy but won't cover majority of array.
        # Possiblity to occur on most distal end.
        # |9997000020|
        valid_length = 0
        for st, end, score in df_dips.select("st", "end", "score").iter_rows():
            if itree_rm[chrom].overlaps(st, end) and score == 0.0:
                valid_length += end - st

        if total_length > 0 and ((valid_length / total_length) > thr_valid_prop):
            chrom, coords = chrom.rsplit(":", 1)
            st, end = coords.split("-")
            # Trim to repeats. Should be slopped afterwards.
            if args.trim_to_repeats:
                _, _, trim_st, trim_end, _ = (
                    group_by_dst(
                        df_rm.filter(pl.col("rtype").is_in(args.trim_to_repeats)),
                        args.dst,
                        "group",
                    )
                    .drop_nulls()
                    .group_by(["group"])
                    .agg(
                        pl.col("chrom").first(), pl.col("st").min(), pl.col("end").max()
                    )
                    .with_columns(len=pl.col("end") - pl.col("st"))
                    .filter(pl.col("len") == pl.col("len").max())
                    .row(0)
                )
                st = int(st)
                new_st = trim_st + st
                new_end = trim_end + st
                print(f"{chrom}\t{new_st}\t{new_end}\t{chrom}:{st}-{end}")
            else:
                print(f"{chrom}\t{st}\t{end}\t{chrom}:{st}-{end}")


if __name__ == "__main__":
    raise SystemExit(main())
