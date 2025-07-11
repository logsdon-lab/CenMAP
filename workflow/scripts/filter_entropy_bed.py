import sys
import argparse
import polars as pl
import scipy.signal
import intervaltree as it
from collections import defaultdict


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
        default=0.5,
        help="Proportion of dips in entropy over all ALR containing regions required to be valid. Majority by default.",
    )
    args = ap.parse_args()

    rm = args.repeatmasker
    # rm = "/project/logsdon_shared/projects/GIAB_HG008_TN/results/6-repeatmasker/repeats/HG008-N/HG008-N_rc-chr14_haplotype2-0000120:4461582-6946730.fa.out"
    bed = args.infile
    # bed = "/project/logsdon_shared/projects/GIAB_HG008_TN/exp_plot/HG008-N_rc-chr14_haplotype2-0000120:4461582-6946730.bed"
    thr_valid_prop = args.prop_valid

    df_rm = pl.read_csv(
        rm,
        separator="\t",
        has_header=False,
        columns=[4, 5, 6, 9],
        new_columns=["chrom", "st", "end", "rtype"],
    ).filter(pl.col("rtype") == "ALR/Alpha")

    itree_rm = defaultdict(it.IntervalTree)
    for chrom, st, end, _ in df_rm.iter_rows():
        itree_rm[chrom].add(it.Interval(st, end))

    df_entropy = (
        pl.read_csv(
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
        .sort(by=["chrom", "st"])
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

        valid_prop = valid_length / total_length
        if valid_prop > thr_valid_prop:
            chrom, coords = chrom.split(":")
            st, end = coords.split("-")
            print(f"{chrom}\t{st}\t{end}")


if __name__ == "__main__":
    raise SystemExit(main())
