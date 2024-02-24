import sys
import argparse
import polars as pl
import edit_distance
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
RGX_CHR = "(chr[0-9XY]+)"


def jaccard_index(a: set[str], b: set[str]) -> float:
    """
    Jaccard similarity index.
    * https://www.statisticshowto.com/jaccard-index/
    """
    return (len(a.intersection(b)) / len(a.union(b))) * 100.0


def format_rm_output(input_path: str) -> pl.LazyFrame:
    return (
        pl.scan_csv(input_path, separator="\t", has_header=False, new_columns=RM_COLS)
        .filter(
            (pl.col("rClass") == "Satellite/centr") | (pl.col("rClass") == "Satellite")
        )
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
        # .with_columns(
        #     start2=pl.when(pl.col("C")=="C").then(pl.col("end")).otherwise(pl.col("start")),
        #     end2=pl.when(pl.col("C")=="C").then(pl.col("start")).otherwise(pl.col("end")),
        # )
        .drop("div", "deldiv", "insdiv", "x", "y", "z", "left", "right", "idx")
    )


def check_cens_status(
    input_rm: str,
    output: TextIO,
    reference_rm: str,
    *,
    match_perc_thr: float = 0.7,
    dst_perc_thr: float = 0.3,
) -> int:
    df_ctg = format_rm_output(input_rm).collect()
    df_ref = (
        format_rm_output(reference_rm)
        .filter(~pl.col("contig").str.starts_with("chm1"))
        .collect()
    )

    contigs, refs, dsts, matches, orts = [], [], [], [], []
    jcontigs, jrefs, jindex = [], [], []
    for ctg, df_ctg_grp in df_ctg.group_by(["contig"]):
        ctg = ctg[0]
        for ref, df_ref_grp in df_ref.group_by(["contig"]):
            ref = ref[0]
            dst_fwd, mtch_fwd = edit_distance.edit_distance(
                df_ref_grp["type"].to_list(), df_ctg_grp["type"].to_list()
            )
            dst_rev, mtch_rev = edit_distance.edit_distance(
                df_ref_grp["type"].to_list(), df_ctg_grp["type"].reverse().to_list()
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
            matches.append(mtch_fwd)
            matches.append(mtch_rev)

    jindex_res = (
        pl.LazyFrame({"contig": jcontigs, "ref": jrefs, "similarity": jindex})
        .with_columns(
            pl.col("similarity").max().over("contig").alias("highest_similarity")
        )
        .filter(pl.col("similarity") == pl.col("highest_similarity"))
        .select(["contig", "ref", "similarity"])
        .collect()
    )
    edit_distance_res = (
        pl.LazyFrame(
            {"contig": contigs, "ref": refs, "dst": dsts, "mtch": matches, "ort": orts}
        )
        # https://stackoverflow.com/a/74630403
        # Calculate percentiles
        .with_columns(
            dst_perc=(pl.col("dst").rank() / pl.col("dst").count()).over("contig"),
            mtch_perc=(pl.col("mtch").rank() / pl.col("mtch").count()).over("contig"),
        )
        # Filter results so only:
        # * Matches gt x percentile.
        # * Distances lt y percentile.
        .filter(
            (pl.col("mtch_perc") > match_perc_thr) & (pl.col("dst_perc") < dst_perc_thr)
        )
        # https://stackoverflow.com/a/74336952
        .with_columns(pl.col("dst").min().over("contig").alias("lowest_dst"))
        .filter(pl.col("dst") == pl.col("lowest_dst"))
        .select(["contig", "ref", "dst", "mtch", "ort"])
        .collect()
    )
    (
        # Join result df
        jindex_res.join(edit_distance_res, on="contig")
        .select(
            contig=pl.col("contig"),
            # Extract chromosome name.
            final_chr=pl.when(pl.col("ref") == pl.col("ref_right"))
            .then(pl.col("ref").str.extract(RGX_CHR))
            # TODO: Fix later. Weird. Cannot use literal string since converted to colname.
            .otherwise(0),
            reorient=pl.col("ort"),
        )
        # Replace chr name in original contig.
        .with_columns(
            final_chr=pl.when(pl.col("final_chr") != "0")
            .then(pl.col("contig").str.replace(RGX_CHR, pl.col("final_chr")))
            .otherwise(pl.col("contig"))
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
        "--match_perc_thr",
        default=0.7,
        type=float,
        help="Matches edit distance percentile threshold. Higher is more stringent.",
    )
    ap.add_argument(
        "--dst_perc_thr",
        default=0.3,
        type=float,
        help="Edit distance percentile threshold. Lower is more stringent.",
    )
    args = ap.parse_args()

    return check_cens_status(
        args.input,
        args.output,
        args.reference,
        match_perc_thr=args.match_perc_thr,
        dst_perc_thr=args.dst_perc_thr,
    )


if __name__ == "__main__":
    raise SystemExit(main())
