import sys
import json
import argparse
from typing import Any
import polars as pl

CT_HEX = "#DDDDDD"
DEF_PATTERNS = {
    "asat": {"pattern": r"ALR", "color": "#58245B"},
    "bsat": {"pattern": r"BSR", "color": "#3A3A3A"},
    "gsat": {"pattern": r"GSAT", "color": "#3A3A3A"},
    "hsat1A": {"pattern": r"SAR", "color": "#3A3A3A"},
    "hsat2": {"pattern": r"^HSATII$", "color": "#3A3A3A"},
    "hsat1B": {"pattern": r"^HSATI$", "color": "#3A3A3A"},
    "hsat3": {"pattern": r"(CATTC)n|(GAATG)n", "color": "#3A3A3A"},
    "ct": {"pattern": None, "color": CT_HEX},
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
RGX_COORDS = r"^(?<ctg>.+):(?<ctg_st>\d+)-(?<ctg_end>\d+)$"


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
        help=(
            "JSON file of regex patterns to label satellite types. "
            "Any type not in patterns is removed. "
            "Requires format: {'satellite': {'pattern': ..., 'color': ...}}"
        ),
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
    ap.add_argument(
        "--add_ct",
        action="store_true",
        help=f"Add CT track with color {CT_HEX} filling entire region. Region corresponds to regex pattern: ({RGX_COORDS})",
    )
    ap.add_argument("--to_abs", action="store_true", help="Convert to absolute coordinates.")

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
    if args.to_abs:
        df = (
            df.with_columns(
                ctg_st=pl.col("chrom").str.extract(r":(\d+)-").fill_null(0).cast(pl.Int64),
            )
            .with_columns(
                chrom_st=pl.col("chrom_st") + pl.col("ctg_st"),
                chrom_end=pl.col("chrom_end") + pl.col("ctg_st"),
            )
        )
    df = (
        df.with_columns(
            strand=pl.when(pl.col("strand") == "C")
            .then(pl.lit("-"))
            .otherwise(pl.lit("+")),
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

    # Remove repeats not matched.
    # Add ct to fill in the rest.
    expr_name = pl.when(pl.lit(False)).then(None)
    expr_item_rgb = pl.when(pl.lit(False)).then(None)
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
    expr_name = expr_name.otherwise(pl.col("rtype"))
    expr_item_rgb = expr_item_rgb.otherwise(None)

    df = df.with_columns(name=expr_name, item_rgb=expr_item_rgb).select(OUTPUT_COLS)

    if args.add_ct:
        df = df.drop_nulls()
        df_ct = (
            df.get_column("chrom")
            .unique()
            .str.extract_groups(RGX_COORDS)
            .to_frame(name="mtch_ctg")
            .unnest("mtch_ctg")
            .with_columns(
                chrom=pl.col("ctg") + ":" + pl.col("ctg_st") + "-" + pl.col("ctg_end"),
            )
            .cast({"ctg_st": pl.Int64, "ctg_end": pl.Int64})
            .with_columns(
                chrom_st=pl.col("ctg_st"),
                chrom_end=pl.col("ctg_end"),
                name=pl.lit("ct"),
                score=pl.lit(0),
                strand=pl.lit("."),
                thick_st=pl.col("ctg_st"),
                thick_end=pl.col("ctg_end"),
                item_rgb=pl.lit(CT_HEX),
            )
            .select(OUTPUT_COLS)
        )
        df = pl.concat((df, df_ct)).sort(by=["chrom", "chrom_st"])
    else:
        df = df.with_columns(pl.col("item_rgb").fill_null(pl.lit(CT_HEX)))

    df.write_csv(args.output_bed, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())
