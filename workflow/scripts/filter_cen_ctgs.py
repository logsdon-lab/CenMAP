import sys
import argparse
import pandas as pd
from io import IOBase
from typing import Any, Iterable, TextIO, TYPE_CHECKING

if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any

DEF_DNA_BRNN_COLS = ("name", "repeat_start", "repeat_end", "repeat_type")
DEF_DNA_BRNN_NAME_SPLIT_COLS = ("ctg_label", "ctg_num", "ctg_start", "ctg_stop")

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
    return None


def add_filt_dna_brnn_args(sub_ap: SubArgumentParser) -> None:
    ap = sub_ap.add_parser(
        "filtdnabrnn",
        description="Script to automatically filter valid repeat contigs from dna-brnn's output.",
    )
    ap.add_argument(
        "-i",
        "--input",
        default=sys.stdin,
        type=argparse.FileType("r"),
        help='Input bed file. Expected to have ("name", "start", "end", "type") columns.',
    )
    ap.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        help="Output bed file.",
        default=sys.stdout,
    )
    ap.add_argument(
        "-c", "--chr", required=True, type=str, help="Chromosome to filter for."
    )
    ap.add_argument(
        "--forward",
        action="store_true",
        help="If sequence orientation is forward (+). Alters repeat length calculation.",
    )
    ap.add_argument(
        "-ci",
        "--columns_in",
        nargs="+",
        help="Input bed cols for dna-nn output. Defaults to (name, repeat_start, repeat_end, repeat_type)",
        default=DEF_DNA_BRNN_COLS,
    )
    ap.add_argument(
        "-cs",
        "--columns_split",
        nargs="+",
        help="Cols created on splitting `name` field by `-` and `:`. Defaults to (ctg_label, ctg_num, ctg_start, ctg_end)",
        default=DEF_DNA_BRNN_NAME_SPLIT_COLS,
    )
    ap.add_argument(
        "--repeat_type",
        type=int,
        default=2,
        help="Repeat type to filter. Defaults to 2 for alpha satellite.",
    )
    ap.add_argument(
        "--repeat_gt_length",
        type=int,
        default=1000,
        help="Repeat length filter. Must be greater than value.",
    )
    return None


def first_line_sv(fh: TextIO, *, delim: str = "\t") -> str:
    first_line = next(fh)
    return first_line.split(delim)


def is_header(first_line: str) -> bool:
    # All words in first line must be alphanumeric to be header.
    return first_line.replace("#", "").replace("_", "").isalpha()


def read_bed_df(input: TextIO, *, input_cols: Iterable[str]) -> pd.DataFrame:
    first_line = first_line_sv(input)
    num_cols = len(first_line)
    has_header = is_header(first_line="".join(first_line))

    if has_header:
        # Let pandas infer.
        return pd.read_csv(input.name, sep="\t")
    else:
        num_input_cols = len(input_cols)
        assert (
            num_cols == len(input_cols)
        ), f"Number of cols not equal to input columns. ({num_cols} != {num_input_cols})"
        return pd.read_csv(input.name, sep="\t", header=0, names=input_cols)


def bedminmax(
    input: TextIO | Iterable[TextIO] | pd.DataFrame,
    output_path: str | None,
    *,
    input_cols: Iterable[str] = DEF_BEDMINMAX_IN_COLS,
    output_cols: Iterable[str] = DEF_BEDMINMAX_OUT_COLS,
    grpby_cols: Iterable[str] = DEF_BEDMINMAX_GRP_COLS,
    sortby_cols: Iterable[str] = DEF_BEDMINMAX_SORT_COLS,
) -> pd.DataFrame | None:
    """
    Collapse records in a bed file by grouping by a series
    of provided columns and merging them based on the
    minimum and maximum values in a `start` and `end` columns, respectively.

    ### Args:
    * `input`:
            * Input bed file, bed files, or `pd.DataFrame`.
    * `output_path`:
            * Output bed file or `pd.DataFrame` if `None`.
    * `input_cols`:
            * Input columns.
            * Expects that a `start` and `end` column exist.
    * `output_cols`:
            * Output columns
    * `grpby_cols`:
            * Group by columns.
    * `sortby_cols`:
            * Sort by columns before writing output.

    ### Raises
    * `AssertionError` if `start` or `end` not a column.

    ### Returns:
    * `None` or `pd.DataFrame`
    """
    if isinstance(input, pd.DataFrame):
        bed_df = input
    elif isinstance(input, IOBase):
        bed_df = read_bed_df(input, input_cols=input_cols)
    else:
        bed_df = pd.concat(
            (read_bed_df(i, input_cols=input_cols) for i in input), axis=0
        )

    assert (
        "start" in bed_df.columns and "end" in bed_df.columns
    ), "Missing required start/end cols."
    bed_out = pd.merge(
        bed_df.groupby(grpby_cols).min()["start"],
        bed_df.groupby(grpby_cols).max()["end"],
        left_index=True,
        right_index=True,
    ).reset_index()

    bed_out = bed_out[output_cols].sort_values(by=sortby_cols)
    if output_path:
        bed_out.to_csv(output_path, sep="\t", header=False, index=False)
        return None
    else:
        return bed_out


def filtdnabrnn(
    input_path: TextIO,
    output_path: TextIO,
    *,
    chr: str,
    forward: bool,
    input_cols: Iterable[str] = DEF_DNA_BRNN_COLS,
    name_split_cols: Iterable[str] = DEF_DNA_BRNN_NAME_SPLIT_COLS,
    repeat_type_filter: int = 2,
    repeat_length_filter: int = 1000,
) -> None:
    """
    Filters the identified repeat output of `dna-nn`.

    ### Args
    * `input`:
            * Bed file output from `dna-nn`,
    * `output_path`:
            * Output filtered bed file
    * `chr`:
            * Chromosome to filter
    * `forward`:
            * If sequences are forward oriented.
    * `input_cols`
            * Columns of input file. Expects the following:
                * `name`
                * `repeat_start`
                * `repeat_stop`
                * `repeat_type`
    * `name_split_cols`:
            * Columns created after splitting `name` col by `-` and `:`.
            * Defaults to (`ctg_label`, `ctg_num`, `ctg_start`, `ctg_stop`)
    * `repeat_type_filter`
            * A label 1 on the 4th column indicates the interval is a region of (AATTC)n ;
            label 2 indicates a region of alpha satellites. [1]
            * Defaults to alpha satellites, a value of 2.
    * `repeat_length_filter`
            * Keep repeats this many bases.
            * Defaults to 1000 bases.

    ### Return
    * `None`

    ### Raise
    * `AssertionError` if not a valid `repeat_type_filter`.

    ### Sources:
    * [1] - https://github.com/lh3/dna-nn#applying-a-trained-model
    """
    assert (
        repeat_type_filter == 1 or repeat_type_filter == 2
    ), f"Invalid repeat type filter. {repeat_type_filter}"

    df_bed = read_bed_df(input_path, input_cols=input_cols)
    # Only look at single chr and only take a specific repeat type.
    df_bed = df_bed.loc[
        df_bed["name"].str.contains(chr) & (df_bed["repeat_type"] == repeat_type_filter)
    ]

    # Split name col. EXPECTS 4 COLS BY DEFAULT.
    # |HG00171_chr16_haplotype1-0000003:4-8430174| -> |HG00171_chr16_haplotype1|0000003|4|8430174|
    df_bed[name_split_cols] = df_bed["name"].str.split(":|-", expand=True, regex=True)
    df_bed.drop(columns=["name"], inplace=True)
    df_bed[["ctg_start", "ctg_stop"]] = df_bed[["ctg_start", "ctg_stop"]].astype(
        "int64"
    )

    # |ooo|x|o|
    # fwd: $3+$5, $3+$6
    # rev: $4-$6, $4-$5
    if forward:
        calc_dst_cols = (
            df_bed["ctg_start"] + df_bed["repeat_start"],
            df_bed["ctg_start"] + df_bed["repeat_stop"],
        )
    else:
        calc_dst_cols = (
            df_bed["ctg_stop"] - df_bed["repeat_stop"],
            df_bed["ctg_stop"] - df_bed["repeat_start"],
        )

    cols = {
        "ctg_label": df_bed["ctg_label"] + "-" + df_bed["ctg_num"]
        if "ctg_num" in df_bed.columns
        else df_bed["ctg_label"],
        **dict(zip(("dst_to_repeat", "dst_to_end_repeat"), calc_dst_cols)),
        "repeat_type": df_bed["repeat_type"],
        "repeat_length": df_bed["repeat_stop"] - df_bed["repeat_start"],
    }

    df_bed = pd.DataFrame(cols)

    # Keep desired repeats above a len.
    df_bed = df_bed.loc[df_bed["repeat_length"] > repeat_length_filter]
    df_bed.to_csv(output_path, sep="\t", header=False, index=False)
    return None


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Script to automatically filter centromeric contigs."
    )
    sub_ap = ap.add_subparsers(dest="cmd")
    add_bedminmax_args(sub_ap)
    add_filt_dna_brnn_args(sub_ap)

    args = ap.parse_args()
    if args.cmd == "bedminmax":
        bedminmax(
            args.input,
            args.output,
            input_cols=args.columns_in,
            output_cols=args.columns_out,
            grpby_cols=args.groupby,
            sortby_cols=args.sortby,
        )
    elif args.cmd == "filtdnabrnn":
        filtdnabrnn(
            args.input,
            args.output,
            chr=args.chr,
            forward=args.forward,
            input_cols=args.columns_in,
            name_split_cols=args.columns_split,
            repeat_type_filter=args.repeat_type,
            repeat_length_filter=args.repeat_gt_length,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
