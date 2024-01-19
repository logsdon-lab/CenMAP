#!/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="")

parser.add_argument(
    "-i", "--bed", help="input 4 column bed file", type=str, required=True
)
parser.add_argument(
    "-o",
    "--output",
    help="output file default /dev/stdout",
    type=str,
    default="/dev/stdout",
)
parser.add_argument(
    "-g",
    "--groupby",
    nargs="+",
    help="Group by columns.",
    default=("chr", "length", "name", "orientation"),
)
parser.add_argument(
    "-ci",
    "--columns_in",
    nargs="+",
    help="Input bed cols. Defaults to chr, start, end, length, name, orientation",
    default=("chr", "start", "end", "length", "name", "orientation"),
)
parser.add_argument(
    "-co",
    "--columns_out",
    nargs="+",
    help="Output bed cols. Defaults to chr, start, end, length, name, orientation",
    default=("chr", "start", "end", "length", "name", "orientation"),
)
parser.add_argument(
    "-s",
    "--sortby",
    nargs="+",
    help="Sort by columns.",
    default=("chr", "start"),
)

args = parser.parse_args()

with open(args.bed, mode="r") as bed_fh:
    first_line = next(bed_fh)
    num_cols = len(first_line.split("\t"))
    # All words in first line must be alphanumeric to be header.
    has_header = (
        first_line.replace("\t", "").replace("#", "").replace("_", "").isalpha()
    )

if has_header:
    # Let pandas infer.
    bed_df = pd.read_csv(args.bed, sep="\t")
else:
    num_input_cols = len(args.columns_in)
    assert num_cols == len(
        args.columns_in
    ), f"Number of cols not equal to input columns. ({num_cols} != {num_input_cols})"
    bed_df = pd.read_csv(args.bed, sep="\t", header=0, names=args.columns_in)

bed_out = pd.merge(
    bed_df.groupby(list(args.groupby)).min()["start"],
    bed_df.groupby(list(args.groupby)).max()["end"],
    left_index=True,
    right_index=True,
).reset_index()

bed_out[list(args.columns_out)].sort_values(by=list(args.sortby)).to_csv(
    args.output, sep="\t", header=False, index=False
)
