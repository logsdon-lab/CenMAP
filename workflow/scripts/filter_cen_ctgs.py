import sys
import argparse
import polars as pl
from io import IOBase
from typing import Any, Iterable, TextIO, TYPE_CHECKING

if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any

DEF_BEDMINMAX_IN_COLS = ("chr", "start", "end", "length", "name", "orientation")
DEF_BEDMINMAX_OUT_COLS = ("chr", "start", "end", "length", "name", "orientation")
DEF_BEDMINMAX_GRP_COLS = ("chr", "length", "name", "orientation")
DEF_BEDMINMAX_SORT_COLS = ("chr", "start")


def add_bedminmax_args(sub_ap: SubArgumentParser) -> None:
    ap = sub_ap.add_parser("bedminmax", description="")
    ap.add_argument(
        "-i",
        "--input",
        help="Input bed file. Defaults to stdin.",
        nargs="+",
        default=sys.stdin,
        type=argparse.FileType("r"),
        required=True,
    )
    ap.add_argument(
        "-o",
        "--output",
        help="Output bed file. Defaults to stdout",
        type=argparse.FileType("w"),
        default=sys.stdout,
    )
    ap.add_argument(
        "-g",
        "--groupby",
        nargs="+",
        help="Group by columns.",
        default=DEF_BEDMINMAX_GRP_COLS,
    )
    ap.add_argument(
        "-ci",
        "--columns_in",
        nargs="+",
        help="Input bed cols. Defaults to chr, start, end, length, name, orientation",
        default=DEF_BEDMINMAX_IN_COLS,
    )
    ap.add_argument(
        "-co",
        "--columns_out",
        nargs="+",
        help="Output bed cols. Defaults to chr, start, end, length, name, orientation",
        default=DEF_BEDMINMAX_OUT_COLS,
    )
    ap.add_argument(
        "-s",
        "--sortby",
        nargs="+",
        help="Sort by columns.",
        default=DEF_BEDMINMAX_SORT_COLS,
    )
    ap.add_argument(
        "--start_col",
        help="Start column to get min.",
        default="start",
    )
    ap.add_argument(
        "--end_col",
        help="End column to get max.",
        default="end",
    )
    ap.add_argument(
        "--allow_empty", action="store_true", help="Allow an output df to be empty."
    )
    return None


def first_line_sv(fh: TextIO, *, delim: str = "\t") -> list[str]:
    first_line = next(fh).strip()
    return first_line.split(delim)


def is_header(first_line: str) -> bool:
    # All words in first line must be alphanumeric to be header.
    return first_line.replace("#", "").replace("_", "").isalpha()


def read_bed_df(input: TextIO, *, input_cols: Iterable[str]) -> pl.DataFrame:
    # If empty, return empty df. Jank ik.
    try:
        first_line = first_line_sv(input)
    except StopIteration:
        return pl.DataFrame()

    num_cols = len(first_line)
    has_header = is_header(first_line="".join(first_line))

    if has_header:
        # Let pandas infer.
        return pl.read_csv(input.name, separator="\t")
    else:
        num_input_cols = len(input_cols)
        assert (
            num_cols == len(input_cols)
        ), f"Number of cols not equal to input columns. ({num_cols} != {num_input_cols})"

        return pl.read_csv(
            input.name, separator="\t", has_header=False, new_columns=input_cols
        )


def bedminmax(
    input: TextIO | Iterable[TextIO] | pl.DataFrame,
    output_path: str | None,
    *,
    start_col: str = "start",
    end_col: str = "end",
    input_cols: Iterable[str] = DEF_BEDMINMAX_IN_COLS,
    output_cols: Iterable[str] = DEF_BEDMINMAX_OUT_COLS,
    grpby_cols: Iterable[str] = DEF_BEDMINMAX_GRP_COLS,
    sortby_cols: Iterable[str] = DEF_BEDMINMAX_SORT_COLS,
    allow_empty: bool = False,
) -> pl.DataFrame | None:
    """
    Collapse records in a bed file by grouping by a series
    of provided columns and merging them based on the
    minimum and maximum values in a `start` and `end` columns, respectively.

    ### Args:
    * `input`:
            * Input bed file, bed files, or `pd.DataFrame`.
    * `output_path`:
            * Output bed file or `pd.DataFrame` if `None`.
    * `start_col`:
            * Start column to get min.
    * `end_col`:
            * End column to get max.
    * `input_cols`:
            * Input columns.
    * `output_cols`:
            * Output columns
    * `grpby_cols`:
            * Group by columns.
    * `sortby_cols`:
            * Sort by columns before writing output.
    * `allow_empty`:
            * Allow empty dataframe to be output.
            * Useful if aggregating and collapsing multiple files.

    ### Raises
    * `AssertionError` if `start` or `end` not a column.

    ### Returns:
    * `None` or `pd.DataFrame`
    """
    if isinstance(input, pl.DataFrame):
        bed_df = input
    elif isinstance(input, IOBase):
        bed_df = read_bed_df(input, input_cols=input_cols)
    else:
        dfs = []
        for i in input:
            df = read_bed_df(i, input_cols=input_cols)
            if not df.is_empty():
                dfs.append(df)
        bed_df = pl.concat(dfs) if dfs else pl.DataFrame()

    # Allow empty df as output.
    if bed_df.is_empty() and allow_empty:
        bed_df.write_csv(output_path, separator="\t", include_header=False)
        return None

    assert (
        start_col in bed_df.columns and end_col in bed_df.columns
    ), "Missing required start/end cols."

    bed_out = bed_df.group_by(grpby_cols).agg(
        [pl.col(start_col).min(), pl.col(end_col).max()]
    )

    bed_out = bed_out.select(output_cols).sort(by=sortby_cols)
    if output_path:
        bed_out.write_csv(output_path, separator="\t", include_header=False)
        return None
    else:
        return bed_out


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Script to automatically filter centromeric contigs."
    )
    sub_ap = ap.add_subparsers(dest="cmd")
    add_bedminmax_args(sub_ap)

    args = ap.parse_args()
    if args.cmd == "bedminmax":
        bedminmax(
            args.input,
            args.output,
            input_cols=args.columns_in,
            output_cols=args.columns_out,
            grpby_cols=args.groupby,
            sortby_cols=args.sortby,
            allow_empty=args.allow_empty,
            start_col=args.start_col,
            end_col=args.end_col,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
