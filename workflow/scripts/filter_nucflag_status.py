import sys
import argparse
import polars as pl

MTYPES = (
    "misjoin",
    "deletion",
    "insertion",
    "other_repeat",
    "false_dup",
    "collapse",
    "scaffold",
)


def main():
    ap = argparse.ArgumentParser(description="Filter NucFlag status bed.")
    ap.add_argument("infile", type=str, help="Input NucFlag status bed with header.")
    ap.add_argument(
        "-q", "--lt_qv", type=float, default=None, help="Filter if less than QV."
    )
    ap.add_argument(
        "-t",
        "--contains_type",
        nargs="*",
        default=MTYPES,
        help="Filter if contains misassembly types",
    )
    args = ap.parse_args()

    thr_qv = args.lt_qv if args.lt_qv else 0.0
    contains_type = [pl.col(col).eq(pl.lit(0.0)) for col in args.contains_type]
    expr_filter = pl.col("QV").ge(thr_qv)
    if contains_type:
        expr_filter = expr_filter & pl.all_horizontal(contains_type)
    df_status = pl.read_csv(args.infile, separator="\t", has_header=True)
    df_status.filter(expr_filter).write_csv(
        sys.stdout, separator="\t", include_header=True
    )


if __name__ == "__main__":
    raise SystemExit(main())
