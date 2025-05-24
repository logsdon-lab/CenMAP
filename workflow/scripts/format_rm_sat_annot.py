import sys
import json
import argparse
from typing import Any
import polars as pl

DEF_PATTERNS = {
    "asat": {"pattern": r"ALR", "color": "#522758"},
    "bsat": {"pattern": r"BSR", "color": "#3A3A3A"},
    "gsat": {"pattern": r"GSAT", "color": "#DDDDDD"},
    "hsat1A": {"pattern": r"SAR", "color": "#3A3A3A"},
    "hsat2": {"pattern": r"^HSATII$", "color": "#3A3A3A"},
    "hsat1B": {"pattern": r"^HSATI$", "color": "#3A3A3A"},
    "hsat3": {"pattern": r"(CATTC)n|(GAATG)n", "color": "#3A3A3A"},
    "ct": {"pattern": None, "color": "#d1d3d4"},
}

OUTPUT_COLS = [
    "chrom",
    "chrom_st",
    "chrom_end",
    "name",
    "score",
    "strand",
    "thick_st",
    "thick_end",
    "item_rgb",
]


def main():
    ap = argparse.ArgumentParser(
        description="Create satellite annotation bed from repeatmasker annotations."
    )

    ap.add_argument(
        "-i",
        "--input_rm",
        help="Input repeatmasker annotations.",
        required=True,
        type=argparse.FileType("rb"),
    )
    ap.add_argument(
        "-p",
        "--patterns",
        help="JSON file of regex patterns to label satellite types. Requires format: {'satellite': {'pattern': ..., 'color': ...}}",
        required=False,
        default=DEF_PATTERNS,
    )
    ap.add_argument(
        "-o",
        "--output_bed",
        help="Output BED9 file.",
        default=sys.stdout,
        type=argparse.FileType("wt"),
    )

    args = ap.parse_args()

    df = pl.read_csv(
        args.input_rm,
        columns=[0, 4, 5, 6, 8, 9, 10],
        separator="\t",
        has_header=False,
        new_columns=[
            "score",
            "chrom",
            "chrom_st",
            "chrom_end",
            "strand",
            "rtype",
            "rclass",
        ],
        truncate_ragged_lines=True,
    )
    df = (
        df.with_columns(
            ctg_st=pl.col("chrom").str.extract(r":(\d+)-").fill_null(0).cast(pl.Int64),
            strand=pl.when(pl.col("strand") == "C")
            .then(pl.lit("-"))
            .otherwise(pl.lit("+")),
        )
        .with_columns(
            chrom_st=pl.col("chrom_st") + pl.col("ctg_st"),
            chrom_end=pl.col("chrom_end") + pl.col("ctg_st"),
        )
        .with_columns(
            thick_st=pl.col("chrom_st"),
            thick_end=pl.col("chrom_end"),
            score=pl.lit(0),
        )
    )

    if isinstance(args.patterns, str):
        with open(args.patterns, "rt") as fh:
            patterns: dict[str, Any] = json.load(fh)
    else:
        patterns = args.patterns

    expr_name = pl.when(pl.lit(False)).then(pl.lit("ct"))
    expr_item_rgb = pl.when(pl.lit(False)).then(pl.lit("#d1d3d4"))
    for sat, pat in patterns.items():
        rgx_pattern = pat.get("pattern")
        color = pat.get("color")

        if not rgx_pattern or not color:
            continue

        expr_name = expr_name.when(pl.col("rtype").str.contains(rgx_pattern)).then(
            pl.lit(sat)
        )
        expr_item_rgb = expr_item_rgb.when(
            pl.col("rtype").str.contains(rgx_pattern)
        ).then(pl.lit(color))
    expr_name = expr_name.otherwise(pl.lit("ct"))
    expr_item_rgb = expr_item_rgb.otherwise(pl.lit("#d1d3d4"))

    df = df.with_columns(name=expr_name, item_rgb=expr_item_rgb).select(OUTPUT_COLS)

    df.write_csv(args.output_bed, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())
