import sys
import re
import argparse
import polars as pl
import editdistance
from typing import TextIO


RM_COLS = [
    "idx",
    "div",
    "deldiv",
    "insdiv",
    "contig",
    "start",
    "end",
    "left",
    "C",
    "type",
    "rClass",
    "right",
    "x",
    "y",
    "z",
]
RGX_CHR = r"(chr[0-9XY]+)"
EDGE_LEN = 100_000


def jaccard_index(a: set[str], b: set[str]) -> float:
    """
    Jaccard similarity index.
    * https://www.statisticshowto.com/jaccard-index/
    """
    return (len(a.intersection(b)) / len(a.union(b))) * 100.0


def format_rm_output(input_path: str) -> pl.LazyFrame:
    return (
        pl.scan_csv(input_path, separator="\t", has_header=False, new_columns=RM_COLS)
        .with_columns(
            pl.col("type")
            .str.replace_many(
                [
                    "/ERVK",
                    "/ERVL",
                    "/ERV1",
                    "/CR1",
                    "/L1",
                    "/L2",
                    "/RTE-X",
                    "/RTE-BovB",
                    "/Gypsy",
                    "-MaLR",
                    "/Alu",
                    "/Deu",
                    "/MIR",
                    "?",
                    "/hAT",
                    "/hAT-Blackjack",
                    "/hAT-Charlie",
                    "/MULE-MuDR",
                    "/PiggyBac",
                    "/TcMar-Mariner" "/TcMar",
                    "/TcMar?",
                    "/hAT-Tip100" "/TcMar-Tigger" "/Dong-R4" "/tRNA",
                ],
                "",
            )
            .str.replace_many(
                [
                    "DNA-Tc2",
                    "DNA?",
                    "DNA-Blackjack",
                    "DNA-Charlie",
                    "DNA-Tigger",
                    "DNA-Tip100",
                ],
                "DNA",
            )
            .str.replace("GSATX", "GSAT")
            .str.replace("LTR\\S", "LTR")
            .str.replace("SAR", "HSat1A")
            .str.replace("HSAT", "HSat1B")
            .str.replace_many(["HSATII", "(CATTC)n", "(GAATG)n"], "HSat2")
        )
        .with_columns(dst=pl.col("end") - pl.col("start"))
        .drop("div", "deldiv", "insdiv", "x", "y", "z", "left", "right", "idx")
    )


def check_cens_status(
    input_rm: str,
    output: TextIO,
    reference_rm: str,
    *,
    dst_perc_thr: float = 0.3,
    edge_perc_alr_thr: float = 0.7,
    edge_len: int = EDGE_LEN,
) -> int:
    df_ctg = format_rm_output(input_rm).collect()
    df_ref = (
        format_rm_output(reference_rm)
        .filter(~pl.col("contig").str.starts_with("chm1"))
        .collect()
    )

    contigs, refs, dsts, orts = [], [], [], []
    jcontigs, jrefs, jindex = [], [], []
    pcontigs, pstatus = [], []
    for ctg, df_ctg_grp in df_ctg.group_by(["contig"]):
        ctg = ctg[0]
        for ref, df_ref_grp in df_ref.group_by(["contig"]):
            ref = ref[0]
            dst_fwd = editdistance.eval(
                df_ref_grp["type"].to_list(),
                df_ctg_grp["type"].to_list(),
            )
            dst_rev = editdistance.eval(
                df_ref_grp["type"].to_list(),
                df_ctg_grp["type"].reverse().to_list(),
            )

            repeat_type_jindex = jaccard_index(
                set(df_ref_grp["type"]), set(df_ctg_grp["type"])
            )
            jcontigs.append(ctg)
            jrefs.append(ref)
            jindex.append(repeat_type_jindex)

            contigs.append(ctg)
            contigs.append(ctg)
            refs.append(ref)
            refs.append(ref)
            orts.append("fwd")
            orts.append("rev")
            dsts.append(dst_fwd)
            dsts.append(dst_rev)

        # Check if partial centromere based on ALR perc on ends.
        # Check N kbp from start and end of contig.
        ledge = df_ctg_grp.filter(pl.col("end") < df_ctg_grp[0]["end"] + edge_len)
        redge = df_ctg_grp.filter(pl.col("end") > df_ctg_grp[-1]["end"] - edge_len)
        try:
            ledge_perc_alr = (
                ledge.group_by("type")
                .agg(pl.col("dst").sum() / ledge["dst"].sum())
                .filter(pl.col("type") == "ALR/Alpha")
                .row(0)[1]
            )
        except Exception:
            ledge_perc_alr = 0.0
        try:
            redge_perc_alr = (
                redge.group_by("type")
                .agg(pl.col("dst").sum() / redge["dst"].sum())
                .filter(pl.col("type") == "ALR/Alpha")
                .row(0)[1]
            )
        except Exception:
            redge_perc_alr = 0.0

        pcontigs.append(ctg)
        pstatus.append(
            ledge_perc_alr > edge_perc_alr_thr or redge_perc_alr > edge_perc_alr_thr
        )

    df_jindex_res = (
        pl.LazyFrame({"contig": jcontigs, "ref": jrefs, "similarity": jindex})
        .with_columns(
            pl.col("similarity").max().over("contig").alias("highest_similarity")
        )
        .filter(pl.col("similarity") == pl.col("highest_similarity"))
        .select(["contig", "ref", "similarity"])
        .collect()
    )
    df_edit_distance_res = pl.DataFrame(
        {"contig": contigs, "ref": refs, "dst": dsts, "ort": orts}
    ).with_columns(
        dst_perc=(pl.col("dst").rank() / pl.col("dst").count()).over("contig"),
    )

    # Filter results so only:
    # * Matches gt x percentile.
    # * Distances lt y percentile.
    dfs_filtered_edit_distance_res = []
    dfs_filtered_ort_same_chr_res = []

    edit_distance_thr_filter = pl.col("dst_perc") < dst_perc_thr
    edit_distance_highest_dst_filter = pl.col("dst_perc") == pl.col("dst_perc").min()

    rgx_chr = re.compile(RGX_CHR)
    for contig, df_edit_distance_res_grp in df_edit_distance_res.group_by(["contig"]):
        chr_name = re.search(rgx_chr, contig[0]).group()
        # Only look at same chr to determine default ort.
        df_edit_distance_res_same_chr_grp = df_edit_distance_res_grp.filter(
            pl.col("ref").str.contains(f"{chr_name}:")
        )

        df_filter_edit_distance_res_grp = df_edit_distance_res_grp.filter(
            edit_distance_thr_filter
        )
        df_filter_ort_res_same_chr_grp = df_edit_distance_res_same_chr_grp.filter(
            edit_distance_thr_filter
        )

        # If none found, default to highest number of matches.
        if df_filter_edit_distance_res_grp.is_empty():
            df_filter_edit_distance_res_grp = df_edit_distance_res_grp.filter(
                edit_distance_highest_dst_filter
            )

        if df_filter_ort_res_same_chr_grp.is_empty():
            df_filter_ort_res_same_chr_grp = df_edit_distance_res_same_chr_grp.filter(
                edit_distance_highest_dst_filter
            )

        dfs_filtered_edit_distance_res.append(df_filter_edit_distance_res_grp)
        dfs_filtered_ort_same_chr_res.append(df_filter_ort_res_same_chr_grp)

    df_filter_edit_distance_res: pl.DataFrame = pl.concat(
        dfs_filtered_edit_distance_res
    )
    df_filter_ort_same_chr_res: pl.DataFrame = pl.concat(dfs_filtered_ort_same_chr_res)

    df_filter_edit_distance_res = (
        df_filter_edit_distance_res
        # https://stackoverflow.com/a/74336952
        .with_columns(pl.col("dst").min().over("contig").alias("lowest_dst"))
        .filter(pl.col("dst") == pl.col("lowest_dst"))
        .select(["contig", "ref", "dst", "ort"])
    )
    # Get pair with lowest dst to get default ort.
    df_filter_ort_same_chr_res = (
        df_filter_ort_same_chr_res.with_columns(
            pl.col("dst").min().over("contig").alias("lowest_dst")
        )
        .filter(pl.col("dst") == pl.col("lowest_dst"))
        .select(["contig", "ort"])
        .rename({"ort": "ort_same_chr"})
    )
    df_partial_contig_res = pl.DataFrame({"contig": pcontigs, "partial": pstatus})

    (
        # Join result dfs.
        # use partial contig res so always get all contigs.
        df_partial_contig_res.join(
            df_jindex_res.join(df_filter_edit_distance_res, on="contig")
            .group_by("contig")
            .first(),
            on="contig",
            how="left",
        )
        # Add default ort per contig.
        .join(df_filter_ort_same_chr_res, on="contig", how="left")
        .select(
            contig=pl.col("contig"),
            # Extract chromosome name.
            # Both results must concur.
            final_contig=pl.when(pl.col("ref") == pl.col("ref_right"))
            .then(pl.col("ref").str.extract(RGX_CHR))
            .otherwise(pl.col("contig").str.extract(RGX_CHR)),
            # Only use orientation if both agree. Otherwise, replace with best same chr ort.
            reorient=pl.when(pl.col("ref") == pl.col("ref_right"))
            .then(pl.col("ort"))
            .otherwise(None)
            .fill_null(pl.col("ort_same_chr")),
            partial=pl.col("partial"),
        )
        # Replace chr name in original contig.
        .with_columns(
            final_contig=pl.col("contig").str.replace(RGX_CHR, pl.col("final_contig")),
            # Never reorient if reference or chrY (ref doesn't contain chrY)
            reorient=pl.when(
                (pl.col("contig").str.starts_with("chm13"))
                | (pl.col("contig").str.contains("chrY"))
            )
            .then(pl.col("reorient").str.replace("rev", "fwd"))
            .otherwise(pl.col("reorient")),
        )
        .write_csv(output, include_header=False, separator="\t")
    )
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Determines if centromeres are incorrectly oriented/mapped with respect to a reference."
    )
    ap.add_argument(
        "-i",
        "--input",
        help="Input RepeatMasker output. Should contain contig reference. Expects no header.",
        type=str,
        required=True,
    )
    ap.add_argument(
        "-o",
        "--output",
        help="List of contigs with actions required to fix.",
        default=sys.stdout,
        type=argparse.FileType("wt"),
    )
    ap.add_argument(
        "-r",
        "--reference",
        required=True,
        type=str,
        help="Reference RM dataframe.",
    )
    ap.add_argument(
        "--dst_perc_thr",
        default=0.3,
        type=float,
        help="Edit distance percentile threshold. Lower is more stringent.",
    )
    ap.add_argument(
        "--edge_perc_alr_thr",
        default=0.7,
        type=float,
        help="Percent ALR on edges of contig to be considered a partial centromere.",
    )
    ap.add_argument(
        "--edge_len",
        default=EDGE_LEN,
        type=int,
        help="Edge len to calculate edge_perc_alr_thr.",
    )
    args = ap.parse_args()

    return check_cens_status(
        args.input,
        args.output,
        args.reference,
        dst_perc_thr=args.dst_perc_thr,
        edge_len=args.edge_len,
        edge_perc_alr_thr=args.edge_perc_alr_thr,
    )


if __name__ == "__main__":
    raise SystemExit(main())
