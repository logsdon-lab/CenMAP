import sys
import argparse

import polars as pl

from collections import deque
from typing import Any, Callable, Iterable, NamedTuple


ACRO_CHRS = {"chr21", "chr22", "chr13", "chr14", "chr15"}


class Interval(NamedTuple):
    begin: int
    end: int
    data: Any

    def length(self):
        return self.end - self.begin


def cmp(x: Interval, y: Interval) -> bool:
    return True


def merge(x: Interval, y: Interval) -> Interval:
    return Interval(x.begin, y.end, x.data)


def merge_itvs(
    itvs: Iterable[Interval],
    dst: int = 1,
    fn_cmp: Callable[[Interval, Interval], bool] | None = None,
    fn_merge_itv: Callable[[Interval, Interval], Interval] | None = None,
) -> list[Interval]:
    if not fn_cmp:
        fn_cmp = cmp
    if not fn_merge_itv:
        fn_merge_itv = merge

    final_itvs = []
    sorted_itvs = deque(sorted(itvs))
    while sorted_itvs:
        try:
            itv_1 = sorted_itvs.popleft()
        except IndexError:
            break
        try:
            itv_2 = sorted_itvs.popleft()
        except IndexError:
            final_itvs.append(itv_1)
            break
        dst_between = itv_2[0] - itv_1[1]
        passes_cmp = fn_cmp(itv_1, itv_2)
        if dst_between <= dst and passes_cmp:
            sorted_itvs.appendleft(fn_merge_itv(itv_1, itv_2))
        else:
            final_itvs.append(itv_1)
            sorted_itvs.appendleft(itv_2)

    return final_itvs


def fn_cmp(itv_1: Interval, itv_2: Interval) -> bool:
    itv_1_chr = itv_1.data
    itv_2_chr = itv_2.data
    # Allow merge if same chromosome or if acrocentrics.
    return itv_1_chr == itv_2_chr


def main():
    ap = argparse.ArgumentParser("Map cens from alignments.")
    ap.add_argument(
        "-i",
        "--input",
        help="Input alignment bed.",
        default=sys.stdin,
        type=argparse.FileType("rb"),
    )
    ap.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("at"))
    ap.add_argument(
        "--allow_multi_chr_prop",
        type=float,
        default=5.0,
        help="Allow multiple chromosome assignment if accounts for greater than x percent of contig.",
    )
    ap.add_argument(
        "--ignore_multi_chr_assignments",
        nargs="*",
        default=ACRO_CHRS,
        help="Ignore chromosomes when building multiple chromosome assignments",
    )

    args = ap.parse_args()

    df = pl.read_csv(
        args.input, separator="\t", has_header=True, truncate_ragged_lines=True
    )

    for qname, df_query in df.sort(by="query_name").group_by(
        ["query_name"], maintain_order=True
    ):
        qname = qname[0]
        query_length = df_query["query_length"].first()
        qitvs = merge_itvs(
            (
                Interval(
                    itv["query_start"],
                    itv["query_end"],
                    itv["#reference_name"],
                )
                for itv in df_query.select(
                    "query_start", "query_end", "#reference_name"
                ).iter_rows(named=True)
            ),
            dst=1,
            fn_cmp=fn_cmp,
        )

        df_ref_length_prop = (
            pl.DataFrame(
                [(itv.data, itv.length()) for itv in qitvs],
                orient="row",
                schema=["#reference_name", "length"],
            )
            .group_by(["#reference_name"], maintain_order=True)
            .agg(prop=100 * (pl.col("length").sum() / query_length))
            .filter(pl.col("prop") > args.allow_multi_chr_prop)
        )

        # Ignore and just take max if all multi_chr.
        if (
            df_ref_length_prop["#reference_name"]
            .is_in(args.ignore_multi_chr_assignments)
            .all()
        ):
            df_ref_length_prop = df_ref_length_prop.filter(
                pl.col("prop") == pl.col("prop").max()
            )

        chrom = "-".join(
            df_ref_length_prop.get_column("#reference_name")
            .unique(maintain_order=True)
            .to_list()
        )
        if not chrom:
            chrom = (
                df_ref_length_prop.filter(pl.col("prop") == pl.col("prop").max())
                .get_column("#reference_name")
                .first()
            )

        df_strand_max = (
            df_query.group_by(["query_name", "strand"])
            .agg(pl.col("matches").sum())
            .filter(pl.col("matches") == pl.col("matches").max())
            .drop("matches")
        )

        df_query_summary = (
            df_query.select("query_name", "query_length")
            .unique()
            .join(df_strand_max, on="query_name")
            .with_columns(
                # Use SRF determined start and end. Not alignment based.
                new_chrom=pl.lit(chrom),
                is_rc=pl.when(pl.col("strand").mode().first() == pl.lit("+"))
                .then(pl.lit("false"))
                .otherwise(pl.lit("true")),
            )
            .select("query_name", "new_chrom", "is_rc", "query_length")
        )
        df_query_summary.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    main()
