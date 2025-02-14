import sys
import json
import argparse
import polars as pl


REPLACE_PATTERNS = dict(
    [
        ("/ERVK", ""),
        ("/ERVL", ""),
        ("/ERV1", ""),
        ("/CR1", ""),
        ("/L1", ""),
        ("/L2", ""),
        ("/RTE-X", ""),
        ("/RTE-BovB", ""),
        ("/Gypsy", ""),
        ("-MaLR", ""),
        ("/Alu", ""),
        ("/Deu", ""),
        ("/MIR", ""),
        ("?", ""),
        ("/hAT-Blackjack", ""),
        ("/hAT-Charlie", ""),
        ("/hAT-Tip100", ""),
        ("/MULE-MuDR", ""),
        ("/PiggyBac", ""),
        ("/TcMar-Mariner", ""),
        ("/TcMar-Tigger", ""),
        ("/TcMar?", ""),
        ("/Dong-R4", ""),
        ("/tRNA-Deu", ""),
        ("DNA-Tc2", "DNA"),
        ("DNA?", "DNA"),
        ("DNA-Blackjack", "DNA"),
        ("DNA-Charlie", "DNA"),
        ("DNA-Tigger", "DNA"),
        ("DNA-Tip100", "DNA"),
        ("GSATX", "GSAT"),
        ("LTR\\S", "LTR"),
        ("SAR", "HSat1A"),
        ("HSATII", "HSat2"),
        ("HSAT", "HSat1B"),
        ("(CATTC)n", "HSat2"),
        ("(GAATG)n", "HSat2"),
    ]
)

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
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i",
        "--infile",
        required=True,
        type=argparse.FileType("rb"),
        help="Input tab-delimited repeatmasker outfile with no header.",
    )
    ap.add_argument(
        "-o",
        "--outfile",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="Output bed9 file.",
    )
    ap.add_argument(
        "-m",
        "--hex_color_mapping",
        required=True,
        type=argparse.FileType("rb"),
        help="Input color mapping JSON.",
    )
    ap.add_argument("-c", "--chrom", default=None, help="Regex for chrom.")
    args = ap.parse_args()

    df = pl.read_csv(
        args.infile,
        separator="\t",
        columns=[0, 4, 5, 6, 8, 9, 10],
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
    )
    rclass_colors = json.load(args.hex_color_mapping)

    if args.chrom:
        df = df.filter(pl.col("chrom").str.contains(args.chrom, literal=False))

    df = (
        df.with_columns(
            ctg_st=pl.col("chrom").str.extract(r":(\d+)-").fill_null(0).cast(pl.Int64),
            ctg_end=pl.col("chrom").str.extract(r"-(\d+)$").fill_null(0).cast(pl.Int64),
        )
        .with_columns(
            chrom_st=pl.col("chrom_st") + pl.col("ctg_st"),
            chrom_end=pl.col("chrom_end") + pl.col("ctg_st"),
            new_rclass=pl.when(pl.col("rclass").str.contains("Satellite"))
            .then(pl.col("rtype"))
            .otherwise(pl.col("rclass")),
        )
        .with_columns(
            new_rclass=pl.col("new_rclass")
            .str.extract(r"(.*?)/")
            .fill_null(pl.col("new_rclass")),
            thick_st=pl.col("chrom_st"),
            thick_end=pl.col("chrom_end"),
        )
        .with_columns(
            item_rgb=pl.col("new_rclass").replace(rclass_colors, default="#FFFFFF"),
            strand=pl.when(pl.col("strand") == "C")
            .then(pl.lit("-"))
            .otherwise(pl.lit("+")),
        )
        .rename({"new_rclass": "name"})
        .select(OUTPUT_COLS)
    )
    df.write_csv(args.outfile, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())
