import sys
import argparse

import polars as pl

from collections import deque
from intervaltree import Interval
from typing import Callable, Iterable


ACRO_CHRS = {"chr21", "chr22", "chr13", "chr14", "chr15"}


def merge_itvs(
    itvs: Iterable[Interval],
    dst: int = 1,
    fn_cmp: Callable[[Interval, Interval], bool] | None = None,
    fn_merge_itv: Callable[[Interval, Interval], Interval] | None = None,
) -> list[Interval]:
    if not fn_cmp:
        fn_cmp = lambda x, y: True
    if not fn_merge_itv:
        fn_merge_itv = lambda x, y: Interval(x.begin, y.end, x.data)

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
        default=33.0,
        help="Allow multiple chromosome assignment if accounts for greater than x percent of contig.",
    )
    ap.add_argument(
        "--ignore_multi_chr_assignments",
        nargs="*",
        default=ACRO_CHRS,
        help="Ignore chromosomes when building multiple chromosome assignments"
    )

    args = ap.parse_args()

    df = (
        pl.read_csv(
            args.input, separator="\t", has_header=True, truncate_ragged_lines=True
        )
        .with_columns(
            query_contig_start=pl.col("query_name").str.extract(r":(\d+)-").fill_null(0),
            query_contig_end=pl.col("query_name").str.extract(r"-(\d+)$").fill_null(0)
        )
        .cast({"query_contig_start": pl.Int64, "query_contig_end": pl.Int64})
        .with_columns(
            pl.col("query_start") + pl.col("query_contig_start"),
            pl.col("query_end") + pl.col("query_contig_start"),
        )
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
                for itv in df_query.select("query_start", "query_end", "#reference_name").iter_rows(named=True)
            ),
            dst=1,
            fn_cmp=fn_cmp,
        )

        df_ref_length_prop = (
            pl.DataFrame(
                [(itv.data, itv.length()) for itv in qitvs],
                orient="row",
                schema=["#reference_name", "length"]
            )
            .group_by(["#reference_name"], maintain_order=True)
            .agg(prop=100 * (pl.col("length").sum() / query_length))
        )
        chrom = "-".join(
            df_ref_length_prop
            .filter((pl.col("prop") > args.allow_multi_chr_prop) & (~pl.col("#reference_name").is_in(args.ignore_multi_chr_assignments)))
            .get_column("#reference_name")
            .unique(maintain_order=True)
            .to_list()
        )
        if not chrom:
            chrom = df_ref_length_prop.filter(pl.col("prop") == pl.col("prop").max()).get_column("#reference_name").first()
        
        df_query_summary = (
            df_query
            .group_by(["query_name"])
            .agg(
                # Use SRF determined start and end. Not alignment based.
                st=pl.col("query_contig_start").first(),
                end=pl.col("query_contig_end").last(),
                new_chrom=pl.lit(chrom),
                is_rc=pl.when(
                    pl.col("strand").mode().first() == pl.lit("+")
                )
                .then(pl.lit("false"))
                .otherwise(pl.lit("true")),
            )
            .with_columns(
                query_name=pl.col("query_name").str.extract(r"^(.+):").fill_null(pl.col("query_name")),
                length=pl.col("end")-pl.col("st")
            )
            .select(
                "query_name", "st", "end", "new_chrom", "length", "is_rc"
            )
        )
        df_query_summary.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    main()
