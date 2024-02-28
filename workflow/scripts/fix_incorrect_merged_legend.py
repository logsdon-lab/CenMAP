import csv
import sys
import argparse


def main() -> int:
    ap = argparse.ArgumentParser(description="Fix incorrect mapped centromeres.")
    ap.add_argument(
        "-ic",
        "--input_corrections",
        help="Input centromere correction lists.",
        required=True,
        type=str,
    )
    ap.add_argument(
        "-l",
        "--input_merged_legend",
        help="Input legend mapping contig name (haplotype-#) to a more complex id (sm_haplotype-#_...).",
        required=True,
        type=str,
    )
    ap.add_argument(
        "-o",
        "--output_agg_legend",
        help="Output aggregated, corrected legend with names replaced from corrections.",
        default=sys.stdout,
        type=argparse.FileType("wt"),
    )

    args = ap.parse_args()

    cens_renamed = {}
    incomplete_cens = set()

    with open(args.input_corrections) as cens_list_fh:
        reader_cens_renamed = csv.reader(cens_list_fh, delimiter="\t")
        for old, new, ort, is_partial in reader_cens_renamed:
            if ort == "rev":
                old = old.replace("chr", "rc_chr")
                new = new.replace("chr", "rc_chr")

            cens_renamed[old] = new

            if is_partial == "true":
                incomplete_cens.add(old)

    writer_legend = csv.writer(args.output_agg_legend, delimiter="\t")
    with open(args.input_merged_legend) as merged_legend_fh:
        reader_merged_legend = csv.reader(merged_legend_fh, delimiter="\t")
        for old, new in reader_merged_legend:
            # Skip incomplete cens.
            if new in incomplete_cens:
                continue
            new_contig_name = cens_renamed.get(new, new)
            writer_legend.writerow((old, new_contig_name))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
