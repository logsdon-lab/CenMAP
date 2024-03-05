import csv
import sys
import re
import argparse

RGX_CHR = re.compile(r"(chr|cen)([0-9XY]+)")


def main():
    ap = argparse.ArgumentParser("Renamed reference stv bed files wtih key.")
    ap.add_argument("-i", "--input_bed", help="input bed")
    ap.add_argument("-l", "--input_len", help="input lengths key")
    ap.add_argument(
        "-o",
        "--output",
        help="output bed with chr name replaced.",
        default=sys.stdout,
        type=argparse.FileType("wt"),
    )

    args = ap.parse_args()

    with open(args.input_len) as in_len:
        lens = {}
        for line in in_len.readlines():
            chr_name, coord = line.strip().split(":")
            if mtch := re.search(RGX_CHR, chr_name):
                chr_lbl = mtch.group(2)
                lens[f"chr{chr_lbl}"] = coord

        assert len(lens) == 23, "Too few or too many chromosome labels."

    with open(args.input_bed, "rt") as in_bed:
        reader = csv.reader(in_bed, delimiter="\t")
        writer = csv.writer(args.output, delimiter="\t")
        for line in reader:
            chr_name = line[0]
            chr_len = lens[chr_name]
            new_name = f"{chr_name}:{chr_len}"
            line[0] = new_name

            writer.writerow(line)


if __name__ == "__main__":
    raise SystemExit(main())
