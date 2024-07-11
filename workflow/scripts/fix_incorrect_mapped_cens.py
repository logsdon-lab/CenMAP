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
        "-ir", "--input_rm_out", help="Input RepeatMasker output files.", nargs="+"
    )
    ap.add_argument(
        "-o",
        "--output_agg_rm",
        help="Output aggregated, corrected RepeatMasker output.",
        default=sys.stdout,
        type=argparse.FileType("wt"),
    )
    ap.add_argument(
        "--ignore_new_names",
        help="Ignore and don't rename/remap names.",
        action="store_true",
    )

    args = ap.parse_args()

    cens_renamed = {}
    incomplete_cens = set()

    with open(args.input_corrections) as cens_list_fh:
        reader_cens_renamed = csv.reader(cens_list_fh, delimiter="\t")
        for old, new, ort, is_partial in reader_cens_renamed:
            if ort == "rev":
                old = old.replace("chr", "rc-chr")
                new = new.replace("chr", "rc-chr")

            if args.ignore_new_names:
                new = old

            cens_renamed[old] = new
            if is_partial == "true":
                incomplete_cens.add(old)

    writer_rm_out = csv.writer(args.output_agg_rm, delimiter="\t")
    for rm_out in args.input_rm_out:
        with open(rm_out) as rm_out_fh:
            reader_rm_out = csv.reader(rm_out_fh, delimiter="\t")
            for line in reader_rm_out:
                contig_name = line[4]
                # Skip incomplete cens.
                if contig_name in incomplete_cens:
                    continue

                new_contig_name = cens_renamed.get(contig_name, contig_name)
                line[4] = new_contig_name
                writer_rm_out.writerow(line)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
