import argparse
import sys
import polars as pl

from intervaltree import Interval
from typing import Callable, Iterable
from collections import deque


# https://github.com/lh3/srf/blob/e54ca8c8eccf6b1f19428b0f862f2c90575290a0/srfutils.js#L336C8-L336C79
INPUT_SRF_COLS = ("chrom", "st", "en", "motif", "wt_ident", "tlen", "len", "keep")
# https://github.com/lh3/TRF-mod/commit/ce9b54f34d571a3f7047ea46ba647bd424ff81dd
INPUT_MON_COLS = (
    "chrom",
    "motif",
    "st",
    "end",
    "period",
    "copyNum",
    "fracMatch",
    "fracGap",
    "score",
    "entroy",
    "pattern",
)
OUTPUT_COLS = ("chrom", "new_st", "new_en", "len")


def group_by_dst(df: pl.DataFrame, dst: int, group_name: str) -> pl.DataFrame:
    try:
        df = df.drop("index")
    except pl.exceptions.ColumnNotFoundError:
        pass
    return (
        df.with_columns(
            # c1  st1 (end1)
            # c1 (st2) end2
            dst_behind=(pl.col("st") - pl.col("en").shift(1)).fill_null(0),
            dst_ahead=(pl.col("st").shift(-1) - pl.col("en")).fill_null(0),
        )
        .with_row_index()
        .with_columns(
            **{
                # Group HOR units based on distance.
                group_name: pl.when(pl.col("dst_behind").le(dst))
                # We assign 0 if within merge dst.
                .then(pl.lit(0))
                # Otherwise, give unique index.
                .otherwise(pl.col("index") + 1)
                # Then create run-length ID to group on.
                # Contiguous rows within distance will be grouped together.
                .rle_id()
            },
        )
        .with_columns(
            # Adjust groups in scenarios where should join group ahead or behind but given unique group.
            # B:64617 A:52416 G:1
            # B:52416 A:1357  G:2 <- This should be group 3.
            # B:1357  A:1358  G:3
            pl.when(pl.col("dst_behind").le(dst) & pl.col("dst_ahead").le(dst))
            .then(pl.col(group_name))
            .when(pl.col("dst_behind").le(dst))
            .then(pl.col(group_name).shift(1))
            .when(pl.col("dst_ahead").le(dst))
            .then(pl.col(group_name).shift(-1))
            .otherwise(pl.col(group_name))
        )
    )


def no_cmp(x: Interval, y: Interval):
    return True


def merge_itv(x: Interval, y: Interval):
    return Interval(x.begin, y.end, x.data)


def merge_itvs(
    itvs: Iterable[Interval],
    dst: int = 1,
    fn_cmp: Callable[[Interval, Interval], bool] | None = None,
    fn_merge_itv: Callable[[Interval, Interval], Interval] | None = None,
) -> list[Interval]:
    if not fn_cmp:
        fn_cmp = no_cmp
    if not fn_merge_itv:
        fn_merge_itv = merge_itv

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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i",
        "--input_bed",
        nargs="+",
        type=argparse.FileType("rb"),
        help=f"Input srf bed. Format: {INPUT_SRF_COLS}",
    )
    ap.add_argument(
        "-m",
        "--monomers_tsv",
        nargs="+",
        type=argparse.FileType("rb"),
        help=f"Input trf monomers. Format: {INPUT_MON_COLS}",
    )
    ap.add_argument(
        "-o",
        "--output_bed",
        type=argparse.FileType("at"),
        default=sys.stdout,
        help=f"Output regions. Format: {OUTPUT_COLS}",
    )
    ap.add_argument(
        "--allowed_monomers",
        type=int,
        nargs="+",
        # 170 - ALR, 42 - SAR
        # https://www.sciencedirect.com/science/article/pii/S1084952122001379
        default=[170, 42],
        help="Allowed monomer lengths.",
    )
    ap.add_argument(
        "--required_monomer",
        type=int,
        default=170,
        help="Required monomer lengths. If required, allowed HORs must also only contain this monomer.",
    )
    ap.add_argument(
        "--top_monomer_by_cn",
        type=int,
        default=3,
        help="Take top n monomers by copy number per chrom and motif.",
    )
    ap.add_argument(
        "--thr_perc_monomer_diff",
        type=float,
        default=1.5,
        help="Percent difference to allow motif based on: 100 * ((motif_len % monomer_length) / monomer_length)",
    )
    ap.add_argument(
        "-g", "--group_by", type=int, default=2_500_000, help="Distance used to group."
    )
    ap.add_argument(
        "-l", "--min_length", type=int, default=30_000, help="Minimum length of region."
    )
    args = ap.parse_args()

    df_srf = pl.concat(
        pl.read_csv(bed, has_header=False, separator="\t", new_columns=INPUT_SRF_COLS)
        for bed in args.input_bed
    )
    df_monomers = (
        pl.concat(
            pl.read_csv(
                tsv, has_header=False, separator="\t", new_columns=INPUT_MON_COLS
            )
            for tsv in args.monomers_tsv
        )
        .with_columns(
            pl.col("chrom")
            .str.extract(r".*?_(rc-chr[0-9XY]+|chr[0-9XY]+)_(?<chrom>.*)", 2)
            .fill_null(pl.col("chrom"))
        )
        .select("chrom", "motif", "period", "copyNum")
        # https://stackoverflow.com/a/78441223
        # Take top 3 motifs with the largest copy number.
        .group_by("chrom", "motif")
        .agg(pl.all().top_k_by("copyNum", k=args.top_monomer_by_cn))
        .explode("period", "copyNum")
    )

    df_monomer_lengths = (
        df_srf.select("chrom")
        .unique()
        .with_columns(mon=args.allowed_monomers)
        .explode("mon")
    )
    df_all_monomers_list = (
        df_monomers.with_columns(period_div=pl.col("period") / args.required_monomer)
        # Do not filter here.
        # Calculate abs difference from required monomer period.
        .with_columns(
            perc_diff=(pl.col("period_div") - pl.col("period_div").round()).abs()
            * 100.0
        )
        .group_by(["chrom", "motif"])
        .agg(all_required_monomers=pl.col("perc_diff") <= args.thr_perc_monomer_diff)
        .with_columns(pl.col("all_required_monomers").list.all())
    )

    df_all_motifs = (
        df_srf.select("chrom", "motif")
        .unique()
        .join(df_monomers, on=["chrom", "motif"])
        .join(df_monomer_lengths, on="chrom")
        .filter(pl.col("period") >= pl.col("mon"))
        .with_columns(period_div=pl.col("period") / pl.col("mon"))
        .with_columns(
            perc_diff=(pl.col("period_div") - pl.col("period_div").round()).abs()
            * 100.0
        )
        .with_columns(valid=pl.col("perc_diff") <= args.thr_perc_monomer_diff)
    )
    df_valid_motifs = df_all_motifs.group_by(["chrom", "motif"]).agg(
        valid=pl.col("valid").any(), mons=pl.col("mon")
    )

    df_srf = (
        df_srf.join(df_valid_motifs, how="left", on=["chrom", "motif"])
        .with_columns(pl.col("valid").fill_null(False))
        .filter(pl.col("valid"))
        # Store if interval contains required monomer
        .with_columns(
            contains_required_monomer=pl.col("mons").list.contains(
                args.required_monomer
            )
        )
        .join(df_all_monomers_list, how="left", on=["chrom", "motif"])
        # Then check if contains a required monomer period, and that all monomers for that HOR must have it.
        # Alpha-satellite HORs should not contain other tandem repeat monomer period aside from ~170 bp.
        # But will keep other monomers if provided (ex. 42bp for SAR/HSAT-1A)
        .filter(
            pl.when(pl.col("contains_required_monomer"))
            .then(pl.col("all_required_monomers"))
            .otherwise(True)
        )
    )

    for chrom, df_chrom in df_srf.sort(by="chrom").group_by(
        ["chrom"], maintain_order=True
    ):
        chrom = chrom[0]
        df_chrom_srf = (
            group_by_dst(df_chrom, args.group_by, "group")
            .group_by(["group"])
            .agg(
                pl.col("chrom").first(),
                pl.col("st").min(),
                pl.col("en").max(),
            )
            .with_columns(rpt_len=pl.col("en") - pl.col("st"))
            .drop("group")
            .sort(by="st")
            .filter(pl.col("rpt_len") > args.min_length)
            .with_columns(
                pl.col("chrom")
                .str.splitn(":", 2)
                .struct.rename_fields(["chrom", "coords"])
            )
            .unnest("chrom")
            .with_columns(
                pl.col("coords")
                .str.splitn("-", 2)
                .struct.rename_fields(["chrom_st", "chrom_en"])
            )
            .unnest("coords")
            .cast({"chrom_st": pl.Int64, "chrom_en": pl.Int64})
            .fill_null(0)
            # Adjust coordinates by contig
            .with_columns(
                new_st=pl.col("st") + pl.col("chrom_st"),
                new_en=pl.col("en") + pl.col("chrom_st"),
            )
        )

        # Must also contain required monomer.
        valid_itvs = []
        for row in df_chrom_srf.iter_rows(named=True):
            df_qry = df_srf.filter(
                (pl.col("chrom") == chrom)
                & (pl.col("st") >= row["st"])
                & (pl.col("en") <= row["en"])
                & pl.col("contains_required_monomer")
            )
            if df_qry.is_empty():
                continue

            valid_itvs.append(Interval(row["new_st"], row["new_en"], row["chrom"]))

        final_itvs = merge_itvs(
            valid_itvs,
            dst=args.group_by,
            fn_cmp=lambda x, y: x.data == y.data,
        )
        for itv in final_itvs:
            print(
                "\t".join([itv.data, str(itv.begin), str(itv.end), str(itv.length())])
            )


if __name__ == "__main__":
    raise SystemExit(main())
