import argparse
import pandas as pd
from typing import Any, Iterable, TextIO, TYPE_CHECKING


if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any

DNA_BRNN_COLS = ("name", "repeat_start", "repeat_end", "repeat_type")
DEF_BEDMINMAX_IN_COLS = ("chr", "start", "end", "length", "name", "orientation")
DEF_BEDMINMAX_OUT_COLS = ("chr", "start", "end", "length", "name", "orientation")
DEF_BEDMINMAX_GRP_COLS = ("chr", "length", "name", "orientation")
DEF_BEDMINMAX_SORT_COLS = ("chr", "start")


def add_bedminmax_args(sub_ap: SubArgumentParser) -> None:
    ap = sub_ap.add_parser("bedminmax", description="")
    ap.add_argument("-i", "--input", help="Input bed file", type=str, required=True)
    ap.add_argument(
        "-o",
        "--output",
        help="Output bed file. Defaults to /dev/stdout",
        type=str,
        default="/dev/stdout",
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
        required=True,
        type=str,
        help='Input bed file. Expected to have ("name", "start", "end", "type") columns.',
    )
    ap.add_argument("-o", "--output", required=True, type=str, help="Output bed file.")
    ap.add_argument(
        "-c", "--chr", required=True, type=str, help="Chromosome to filter for."
    )
    return None


def first_line_sv(fh: TextIO, *, delim: str = "\t") -> str:
    first_line = next(fh)
    return first_line.split(delim)


def is_header(first_line: str) -> bool:
    # All words in first line must be alphanumeric to be header.
    return first_line.replace("\t", "").replace("#", "").replace("_", "").isalpha()


def read_bed_df(input: str, *, input_cols: Iterable[str]) -> pd.DataFrame:
    with open(input, mode="r") as bed_fh:
        first_line = first_line_sv(bed_fh)
        num_cols = len(first_line)
        has_header = is_header(first_line_sv)

    if has_header:
        # Let pandas infer.
        return pd.read_csv(input, sep="\t")
    else:
        num_input_cols = len(input_cols)
        assert (
            num_cols == len(input_cols)
        ), f"Number of cols not equal to input columns. ({num_cols} != {num_input_cols})"
        return pd.read_csv(input, sep="\t", header=0, names=input_cols)


def bedminmax(
    input: str | pd.DataFrame,
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
            * Input bed file or `pd.DataFrame`.
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
    else:
        bed_df = read_bed_df(input, input_cols=input_cols)

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


# (Per chr)
# TODO: Where does this come from?
# /net/eichler/vol27/projects/hgsvc/nobackups/analysis/glennis/map_align_t2t_chm13/results/t2t_chm13_v2.0_cens/bed/centromeres/dna-brnn/chm13_cens.trimmed.bed
# grep "chr1:" chm13_cens.trimmed.bed | \
# sed 's/:/\t/g' | \
# sed 's/-/\t/g' | \
# awk -v OFS="\t" '{print $1, $2+$4, $2+$5, $6, $5-$4}' | \
# awk '$4==2' | \
# awk '$5>1000' > chr1_tmp.fwd.bed

# (Per chr + sample)
# grep "chr2_" ${sample}.renamed.fwd.bed | \
# sed 's/:/\t/g' | \
# sed 's/-/\t/g' | \
# awk -v OFS="\t" '{print $1"-"$2, $3+$5, $3+$6, $7, $6-$5}' | \
# awk '$4==2' | \
# awk '$5>1000' >> chr2_tmp.fwd.bed


# TODO: Needs to be done after all samples evaluated.
# /net/eichler/vol28/home/glogsdon/utilities/bedminmax.py \
# -i chr2_tmp.fwd.bed | \
# awk -v OFS="\t" '{print $1, $2, $3, $3-$2}' | \
# awk -v OFS="\t" '{print $1, $2-467987, $3+522450, $3-$2}' | \
# awk '$4>1000000' | \
# awk -v OFS="\t" '$2<0 {$2=0}1' > chr2_contigs.fwd.repeat.bed
def filtdnabrnn(
    input_path: str,
    output_path: str,
    *,
    chr: str,
    input_cols: Iterable[str] = DNA_BRNN_COLS,
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
    * `input_cols`
            * Columns of input file. Expects the following:
                * `name` (sample, chr, haplotype, contig, regstart, regstop)
                * `repeat_start`
                * `repeat_stop`
                * `repeat_type`
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
        df_bed["repeat_type"].str.contains(chr)
        & (df_bed["repeat_type"] == repeat_type_filter)
    ]

    # Split name col
    # |HG00171_chr16_haplotype1-0000003:4-8430174| -> |HG00171_chr16_haplotype1|0000003|4|8430174|
    df_bed[["sample_chr_haplo", "ctg_num", "ctg_start", "ctg_stop"]] = df_bed[
        "name"
    ].str.split(":|-", expand=True, regex=True)
    df_bed.drop(columns=["name"], inplace=True)

    # |ooo|x|o|
    output_cols = [
        "name",
        "dst_to_repeat",
        "dst_to_end_repeat",
        "repeat_type",
        "repeat_length",
    ]
    df_bed = pd.concat(
        [
            df_bed["sample_chr_haplo"] + "-" + df_bed["ctg_num"],
            df_bed["ctg_start"] + df_bed["repeat_start"],
            df_bed["ctg_start"] + df_bed["repeat_stop"],
            df_bed["repeat_type"],
            df_bed["repeat_stop"] - df_bed["repeat_start"],
        ],
        names=output_cols,
        axis=1,
    )

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
            args.columns_in,
            args.column_out,
            args.groupby,
            args.sortby,
        )
    elif args.cmd == "filtdnabrnn":
        filtdnabrnn(args.input, args.output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
