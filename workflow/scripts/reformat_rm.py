import sys
import csv
import argparse
import logging

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
    stream=sys.stderr,
)


def main():
    ap = argparse.ArgumentParser("Reformat RepeatMasker output.")

    ap.add_argument("-i", "--infile", required=True, help="RepeatMasker output file.")
    ap.add_argument(
        "-of", "--original_faidx", required=True, help="Original fasta index file."
    )
    ap.add_argument(
        "-rf", "--renamed_faidx", required=True, help="Renamed fasta index file."
    )
    ap.add_argument(
        "-o",
        "--outfile",
        help="Output file.",
        default=sys.stdout,
        type=argparse.FileType("wt"),
    )

    args = ap.parse_args()

    rename_legend = {}
    with (
        open(args.original_faidx, "rt") as o_fai,
        open(args.renamed_faidx, "rt") as r_fai,
    ):
        o_fai_reader = csv.reader(o_fai, delimiter="\t")
        r_fai_reader = csv.reader(r_fai, delimiter="\t")
        for original_line, renamed_line in zip(o_fai_reader, r_fai_reader):
            original_name = original_line[0]
            renamed_name = renamed_line[0]

            rename_legend[renamed_name] = original_name

    with (
        open(args.infile, "rt") as infile,
    ):
        # Skip first four lines.
        for _ in range(3):
            next(infile)
        rm_writer = csv.writer(args.outfile, delimiter="\t")
        for line in infile.readlines():
            line = line.strip().split()
            renamed_name = line[4]
            if renamed_name in rename_legend:
                line[4] = rename_legend[renamed_name]
            else:
                logging.warning(
                    f"Could not find {renamed_name} in rename_legend. Writing original name."
                )

            rm_writer.writerow(line)


if __name__ == "__main__":
    raise SystemExit(main())
