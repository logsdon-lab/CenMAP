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
        help="RepeatMasker output file with no header.",
    )
    ap.add_argument(
        "-l",
        "--thr_length",
        type=int,
        default=0,
        help="Minimum length of centromere.",
    )
    args = ap.parse_args()

    rm = args.repeatmasker
    # rm = "/project/logsdon_shared/projects/GIAB_HG008_TN/results/6-repeatmasker/repeats/HG008-N/HG008-N_rc-chr14_haplotype2-0000120:4461582-6946730.fa.out"
    bed = args.infile
    # bed = "/project/logsdon_shared/projects/GIAB_HG008_TN/exp_plot/HG008-N_rc-chr14_haplotype2-0000120:4461582-6946730.bed"
    thr_length = args.thr_length

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
    ).with_row_index()
    for chrom, df_chrom_entropy in df_entropy.group_by(["chrom"]):
        chrom: str = chrom[0]
        # Find dips.
        dips, _ = scipy.signal.find_peaks(-df_chrom_entropy["score"].to_numpy())
        df_dips = df_chrom_entropy.filter(pl.col("index").is_in(dips))

        # Is dip, contains ALR, and greater than thr_length
        if any(
            itree_rm[chrom].overlaps(st, end) and score == 0.0 and end - st > thr_length
            for st, end, score in df_dips.select("st", "end", "score").iter_rows()
        ):
            chrom, coords = chrom.split(":")
            st, end = coords.split("-")
            print(f"{chrom}\t{st}\t{end}")


if __name__ == "__main__":
    raise SystemExit(main())
