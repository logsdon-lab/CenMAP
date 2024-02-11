import sys
import json
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser(
        description="Get centromere coord/data from dna-brnn output."
    )
    ap.add_argument(
        "-i", "--input", help="Input bed file.", default=sys.stdin, required=True
    )
    ap.add_argument(
        "-o",
        "--output",
        help="Output json file with start and end positions.",
        default=sys.stdout,
    )
    args = ap.parse_args()

    df = pd.read_csv(
        args.input,
        sep="\t",
        header=None,
        names=["lbl", "start", "end", "repeat"],
    )
    cen_data = {}
    chr_name, str_ctg_coords = df["lbl"].iloc[0].split(":")
    ctg_dst = abs(eval(str_ctg_coords))

    # TODO: Ask if start and end are correct terms.
    if chr_name == "chr16":
        # Get row with max diff in start and end
        rep_lens = df["end"] - df["start"]
        row = df.loc[rep_lens.idxmax()]
        cen_data["start"] = int(ctg_dst - row["start"])
        cen_data["end"] = int(ctg_dst - row["end"])
    else:
        cen_data["start"] = int(df["start"].min())
        cen_data["end"] = int(ctg_dst - df["end"].max())

    try:
        print(json.dumps(cen_data), file=args.output)
    except AttributeError:
        with open(args.output, "wt") as fh:
            json.dump(cen_data, fh)


if __name__ == "__main__":
    raise SystemExit(main())
