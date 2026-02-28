import sys
import gzip
import tarfile
import pathlib
import polars as pl

from typing import Any, Callable
from multiprocessing import get_context
from concurrent.futures import Future, ProcessPoolExecutor

OUTPUT_DIR_TYPES = ("bed", "plot")
BED9_COLS = (
    "chrom",
    "st",
    "end",
    "name",
    "score",
    "strand",
    "tst",
    "tend",
    "item_rgb",
)


def write_bedfile(
    source: pathlib.Path,
    df_rename_key: pl.DataFrame,
    output_dir: pathlib.Path,
    output_fname_suffix: str,
    *,
    kwargs_polars: dict[str, Any] | None = None,
    kwargs_format_bedfile: dict[str, Any] | None = None,
    split_by_sample: bool = False,
    fn_df_format: Callable[[pl.DataFrame], pl.DataFrame] | None = None,
):
    kwarg_polars = kwargs_polars if kwargs_polars else {}
    kwargs_format_bedfile = kwargs_format_bedfile if kwargs_format_bedfile else {}
    df_bed = pl.read_csv(
        source,
        # glob=True,
        separator="\t",
        has_header=False,
        # new_columns=BED9_COLS,
        **kwarg_polars,
    )
    df_bed_adj = format_bedfile(
        df_bed,
        df_rename_key,
        **kwargs_format_bedfile,
        # convert_to_abs=True
    )
    if fn_df_format:
        df_bed_adj = fn_df_format(df_bed_adj)

    if split_by_sample:
        dfs_sm = (
            (str(sm[0]), df)
            for sm, df in df_bed_adj.partition_by(by="sm", as_dict=True).items()
        )
    else:
        dfs_sm = (("all", df) for df in [df_bed_adj])
    for sm, df in dfs_sm:
        assert output_fname_suffix.endswith(".gz"), (
            "Output filename suffix must end in gz."
        )
        with gzip.open(
            output_dir.joinpath("bed", f"{sm}_{output_fname_suffix}"), "wb"
        ) as fh:
            df.write_csv(fh, separator="\t", include_header=False)
    return


def format_bedfile(
    df_bed: pl.DataFrame,
    df_rename_key: pl.DataFrame,
    *,
    convert_to_abs: bool = False,
    convert_to_original: bool = True,
) -> pl.DataFrame:
    df_bed = (
        df_bed.with_columns(
            mtch=pl.col("chrom").str.extract_groups(
                r"^(?<ctg>.+):(?<ctg_st>[0-9]+)-[0-9]+$"
            )
        )
        .unnest("mtch")
        .with_columns(
            chrom=pl.col("ctg").fill_null(pl.col("chrom")),
            ctg_st=pl.col("ctg_st").cast(pl.UInt64).fill_null(pl.lit(0)),
        )
        .drop("ctg")
    )

    if convert_to_abs:
        df_bed = df_bed.with_columns(
            st=pl.col("st") + pl.col("ctg_st"),
            end=pl.col("end") + pl.col("ctg_st"),
        )
        if "tst" in df_bed.columns:
            df_bed = df_bed.with_columns(
                tst=pl.col("tst") + pl.col("ctg_st"),
                tend=pl.col("tend") + pl.col("ctg_st"),
            )

    df_bed = df_bed.drop("ctg_st")
    if convert_to_original:
        df_bed = (
            df_bed.join(df_rename_key, left_on="chrom", right_on="new_ctg", how="left")
            .with_columns(
                st=pl.when(pl.col("chrom").str.contains("rc-chr", literal=True))
                .then(pl.col("ctg_len") - pl.col("end"))
                .otherwise(pl.col("st")),
                end=pl.when(pl.col("chrom").str.contains("rc-chr", literal=True))
                .then(pl.col("ctg_len") - pl.col("st"))
                .otherwise(pl.col("end")),
            )
            .cast({"st": pl.UInt64, "end": pl.UInt64})
        )

        if "tst" in df_bed.columns:
            df_bed = df_bed.with_columns(
                tst=pl.when(pl.col("chrom").str.contains("rc-chr", literal=True))
                .then(pl.col("ctg_len") - pl.col("tend"))
                .otherwise(pl.col("tst")),
                tend=pl.when(pl.col("chrom").str.contains("rc-chr", literal=True))
                .then(pl.col("ctg_len") - pl.col("tst"))
                .otherwise(pl.col("tend")),
            ).cast({"tst": pl.UInt64, "tend": pl.UInt64})

        df_bed = (
            df_bed.with_columns(chrom=pl.col("ctg"))
            .drop_nulls()
            # Keep sm
            .drop("ctg", "ctg_len")
        )
    return df_bed


def cp_rm_outputs(
    input_dir: pathlib.Path, output_dir: pathlib.Path, *, split_by_sample: bool = False
) -> None:
    df_rename_key = get_rename_key(input_dir)
    # RepeatMasker
    rm_dir = input_dir.joinpath("6-repeatmasker")
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/6-repeatmasker/bed/chr_all_og.bed
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/6-repeatmasker/plots/chr_all_cens_og.pdf
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/6-repeatmasker/plots/chr_all_cens_og.png

    write_bedfile(
        source=rm_dir.joinpath("bed", "chr_*.bed"),
        df_rename_key=df_rename_key,
        output_dir=output_dir,
        output_fname_suffix="rm.bed.gz",
        kwargs_polars={"glob": True, "new_columns": BED9_COLS},
        split_by_sample=split_by_sample,
    )
    with tarfile.open(output_dir.joinpath("plot", "all_rm.tar.gz"), "w:gz") as tar:
        for afile in rm_dir.joinpath("plots").glob("chr_*_cens_og.*"):
            tar.add(afile, arcname=afile.name)


def cp_rm_fix_outputs(
    input_dir: pathlib.Path, output_dir: pathlib.Path, *, split_by_sample: bool = False
) -> None:
    df_rename_key = get_rename_key(input_dir)
    # RepeatMasker (remove non-ALR and partial centromeres, if any)
    rm_check_dir = next(input_dir.glob("7.1-fix_cens_w_*"))
    rm_plots = rm_check_dir.joinpath("plots").glob("chr_*_cens.*")
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/7.1-fix_cens_w_repeatmasker/bed/chr_all.bed
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/7.1-fix_cens_w_repeatmasker/plots/chr_all_cens.pdf
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/7.1-fix_cens_w_repeatmasker/plots/chr_all_cens.png

    write_bedfile(
        source=rm_check_dir.joinpath("bed", "chr_*.bed"),
        df_rename_key=df_rename_key,
        output_dir=output_dir,
        output_fname_suffix="rm_fix.bed.gz",
        kwargs_polars={"glob": True, "new_columns": BED9_COLS},
        split_by_sample=split_by_sample,
    )
    with tarfile.open(output_dir.joinpath("plot", "all_rm_fix.tar.gz"), "w:gz") as tar:
        for afile in rm_plots:
            tar.add(afile, arcname=afile.name)


def cp_cdr_finder_outputs(
    input_dir: pathlib.Path, output_dir: pathlib.Path, *, split_by_sample: bool = False
) -> None:
    df_rename_key = get_rename_key(input_dir)
    # CDR-Finder
    cdr_finder_dir = input_dir.joinpath("8-cdr_finder")
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/8-cdr_finder/bed/all_cdrs.bed
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/8-cdr_finder/bed/all_binned_freq.bed

    write_bedfile(
        source=cdr_finder_dir.joinpath("bed", "all_cdrs.bed"),
        df_rename_key=df_rename_key,
        output_dir=output_dir,
        output_fname_suffix="cdrs.bed.gz",
        kwargs_polars={"new_columns": ["chrom", "st", "end"]},
        split_by_sample=split_by_sample,
    )
    write_bedfile(
        source=cdr_finder_dir.joinpath("bed", "all_binned_freq.bed"),
        df_rename_key=df_rename_key,
        output_dir=output_dir,
        output_fname_suffix="cpg_methyl.bedgraph.gz",
        kwargs_polars={
            "columns": list(range(4)),
            "new_columns": ["chrom", "st", "end", "avg_cpg_methyl"],
        },
        split_by_sample=split_by_sample,
    )


def cp_nucflag_outputs(
    input_dir: pathlib.Path, output_dir: pathlib.Path, *, split_by_sample: bool = False
) -> None:
    df_rename_key = get_rename_key(input_dir)
    # NucFlag
    nucflag_dir = input_dir.joinpath("8-nucflag")
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/8-nucflag/HG008-N_v6.3_status.bed
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/8-nucflag/HG008-N_v6.3_misassemblies.bed
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/8-nucflag/HG008-N_v6.3 #
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/8-nucflag/HG008-T_v3.2 #

    # df_status_bed = pl.read_csv(
    #     nucflag_dir.joinpath("*_status.bed"),
    #     glob=True,
    #     separator="\t",
    #     has_header=False,
    #     new_columns=["chrom", "st", "end", "name"]
    # )
    # df_status_bed_adj = format_bedfile(df_status_bed, df_rename_key)
    # with gzip.open(output_dir.joinpath("bed", "all_nucflag_status.bed.gz"), "wb") as fh:
    #     df_status_bed_adj.write_csv(fh, separator="\t", include_header=False)

    write_bedfile(
        source=nucflag_dir.joinpath("*_misassemblies.bed"),
        df_rename_key=df_rename_key,
        output_dir=output_dir,
        output_fname_suffix="nucflag_misassemblies.bed.gz",
        kwargs_polars={"glob": True, "new_columns": ["chrom", "st", "end", "name"]},
        split_by_sample=split_by_sample,
    )


def cp_sat_annot_outputs(
    input_dir: pathlib.Path, output_dir: pathlib.Path, *, split_by_sample: bool = False
) -> None:
    df_rename_key = get_rename_key(input_dir)
    # Satellite annotation
    rm_sat_annot_dir = input_dir.joinpath("8-format_repeatmasker_sat_annot")

    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/8-format_repeatmasker_sat_annot/bed/all_cens.annotation.bed
    write_bedfile(
        source=rm_sat_annot_dir.joinpath("bed", "all_cens.annotation.bed"),
        df_rename_key=df_rename_key,
        output_dir=output_dir,
        output_fname_suffix="rm_sat_annot.bed.gz",
        kwargs_polars={"new_columns": BED9_COLS},
        split_by_sample=split_by_sample,
    )


def cp_humas_mon_stv_outputs(
    input_dir: pathlib.Path, output_dir: pathlib.Path, *, split_by_sample: bool = False
) -> None:
    df_rename_key = get_rename_key(input_dir)
    # Alpha-satellite monomer and HOR structural variants bed
    humas_annot_dir = input_dir.joinpath("8-humas_annot")
    humas_annot_fmt_dir = input_dir.joinpath("10-format_hor_stv")

    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/8-humas_annot/*/final_decomposition.bed # monomers_all.bed
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/10-format_hor_stv/bed/all/stv_all.bed
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/10-format_hor_stv/bed/all/stv_complete.bed

    write_bedfile(
        source=humas_annot_dir.joinpath("*", "final_decomposition.bed"),
        df_rename_key=df_rename_key,
        output_dir=output_dir,
        output_fname_suffix="monomers.bed.gz",
        kwargs_polars={"new_columns": BED9_COLS, "glob": True},
        kwargs_format_bedfile={"convert_to_abs": True},
        split_by_sample=split_by_sample,
    )
    write_bedfile(
        source=humas_annot_fmt_dir.joinpath("bed", "all", "stv_all.bed"),
        df_rename_key=df_rename_key,
        output_dir=output_dir,
        output_fname_suffix="hor_stv.bed.gz",
        kwargs_polars={"new_columns": BED9_COLS},
        split_by_sample=split_by_sample,
    )
    write_bedfile(
        source=humas_annot_fmt_dir.joinpath("bed", "all", "stv_complete.bed"),
        df_rename_key=df_rename_key,
        output_dir=output_dir,
        output_fname_suffix="hor_stv_complete.bed.gz",
        kwargs_polars={"new_columns": BED9_COLS},
        split_by_sample=split_by_sample,
    )


def cp_hor_array_length_outputs(
    input_dir: pathlib.Path, output_dir: pathlib.Path, *, split_by_sample: bool = False
) -> None:
    df_rename_key = get_rename_key(input_dir)
    # Alpha-satellite monomer and HOR structural variants bed
    hor_arr_len_dir = input_dir.joinpath("11-calculate_hor_length")

    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/11-calculate_hor_length/bed/all_AS-HOR_lengths.bed
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/11-calculate_hor_length/bed/all_AS-HOR_strand.bed

    write_bedfile(
        source=hor_arr_len_dir.joinpath("bed", "all_AS-HOR_lengths.bed"),
        df_rename_key=df_rename_key,
        output_dir=output_dir,
        output_fname_suffix="as_hor_lengths.bed.gz",
        kwargs_polars={
            "new_columns": ["chrom", "st", "end", "length", "num_hors", "cov_hor"]
        },
        split_by_sample=split_by_sample,
    )
    write_bedfile(
        source=hor_arr_len_dir.joinpath("bed", "all_AS-HOR_strand.bed"),
        df_rename_key=df_rename_key,
        output_dir=output_dir,
        output_fname_suffix="as_hor_lengths_strand.bed.gz",
        kwargs_polars={
            "new_columns": [
                "chrom",
                "st",
                "end",
                "length",
                "num_hors",
                "cov_hor",
                "strand",
            ]
        },
        split_by_sample=split_by_sample,
    )


def cp_moddotplot_outputs(
    input_dir: pathlib.Path, output_dir: pathlib.Path, *, split_by_sample: bool = False
) -> None:
    df_rename_key = get_rename_key(input_dir)
    # ModDotPlot
    moddotplot_dir = input_dir.joinpath("11-moddotplot")

    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/11-moddotplot/combined/plots #
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/11-moddotplot/original/all/HG008-N_v6.3_chr1_hap1:120587249-128026751/HG008-N_v6.3_chr1_hap1:120587249-128026751.bed
    write_bedfile(
        source=moddotplot_dir.joinpath("original", "all", "*", "*.bed"),
        df_rename_key=df_rename_key,
        output_dir=output_dir,
        output_fname_suffix="moddotplot.bedpe.gz",
        # tst and tend is ref start and ref end. Just need to reorient.
        kwargs_polars={
            "glob": True,
            "new_columns": ["chrom", "st", "end", "rchrom", "tst", "tend", "ident"],
        },
        kwargs_format_bedfile={"convert_to_abs": True},
        fn_df_format=lambda df: df.with_columns(rchrom=pl.col("chrom")),
        split_by_sample=split_by_sample,
    )
    with tarfile.open(
        output_dir.joinpath("plot", "all_moddotplot_new_coords.tar.gz"), "w:gz"
    ) as tar:
        for afile in moddotplot_dir.joinpath("combined", "plots", "all").glob("*.p*"):
            tar.add(afile, arcname=afile.name)


def get_rename_key(input_results_dir: pathlib.Path) -> pl.DataFrame:
    dfs_key = []
    for file_key in input_results_dir.joinpath("5-ident_cen_ctgs", "bed").glob(
        "*_rename_key.tsv"
    ):
        sm = file_key.stem.replace("_rename_key", "")
        df_key = pl.read_csv(
            file_key,
            separator="\t",
            has_header=False,
            new_columns=["ctg", "new_ctg", "ctg_len"],
        ).with_columns(sm=pl.lit(sm))
        dfs_key.append(df_key)

    return pl.concat(dfs_key)


def main():
    input_results_dir = pathlib.Path(sys.argv[1])
    # /project/logsdon_shared/projects/HG008_TN/CenMAP_T-N_v3.2-v6.3/results/final
    output_dir = pathlib.Path(sys.argv[2])

    for dtype in OUTPUT_DIR_TYPES:
        output_dir.joinpath(dtype).mkdir(exist_ok=True, parents=True)

    with ProcessPoolExecutor(max_workers=12, mp_context=get_context("spawn")) as pool:
        futures: list[Future] = []
        for fn in (
            cp_rm_outputs,
            cp_rm_fix_outputs,
            cp_cdr_finder_outputs,
            cp_nucflag_outputs,
            cp_sat_annot_outputs,
            cp_humas_mon_stv_outputs,
            cp_hor_array_length_outputs,
            cp_moddotplot_outputs,
        ):
            futures.append(
                pool.submit(fn, input_results_dir, output_dir, split_by_sample=True)
            )

        for future in futures:
            if future.cancelled():
                print(future.exception(), file=sys.stderr)
            else:
                _ = future.result()


if __name__ == "__main__":
    raise SystemExit(main())
