import argparse
import numpy as np
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from collections import OrderedDict


DEF_COLS = ("chrom", "chrom_st", "chrom_end", "length")
DEF_CHR_ORDER = [f"chr{i}" for i in (*range(1, 23), "X", "Y")]
DEF_CHR_COLORS = dict(
    zip(
        DEF_CHR_ORDER,
        (
            "#403E80",
            "#2C477D",
            "#346B9B",
            "#3C80AA",
            "#4587A2",
            "#589D96",
            "#73ACA4",
            "#87BAAF",
            "#94BE9F",
            "#9FC38C",
            "#A1C27C",
            "#A4C165",
            "#C2C969",
            "#B5A957",
            "#DFD06C",
            "#F2D46C",
            "#E3C765",
            "#E5BA61",
            "#D89E56",
            "#C8824A",
            "#BB6B3E",
            "#AB5C40",
            "#9B4C41",
            "#8C3C42",
        ),
    )
)


def main():
    ap = argparse.ArgumentParser(
        description="Plot cumulative centromere HOR array lengths."
    )
    ap.add_argument(
        "-i",
        "--infile",
        help="Input centromere HOR array lengths.",
        type=argparse.FileType("rb"),
    )
    ap.add_argument(
        "-a",
        "--added_inputs",
        help="Additional centromere HOR array lengths. First column should be be a subset of --chroms",
        nargs="*",
        metavar="{lbl}={path}={color}",
        type=str,
    )
    ap.add_argument(
        "-m",
        "--mode",
        type=str,
        choices=["total", "arr"],
        default="total",
        help="Plotting mode. Either total live array length (total) or live array length (arr).",
    )
    ap.add_argument(
        "-c",
        "--chroms",
        nargs="+",
        default=DEF_CHR_ORDER,
        help="Chromosome names.",
    )
    ap.add_argument(
        "--chrom_colors",
        default=None,
        help="Chromosome colors as TSV file with chrom to color mapping.",
    )
    ap.add_argument("-o", "--output", help="Output plot file.", required=True, type=str)
    args = ap.parse_args()

    # Reverse to prevent matching chr1 with both chr1 and chr11
    rgx_chrom = "|".join([*reversed(sorted(args.chroms)), "-"])
    rgx_name_groups = r"^.*?_(?<chrom_name>(" + rgx_chrom + r")*)_.*?$"
    plot_all = "all" in args.chroms

    if args.chrom_colors:
        _chroms = set(args.chroms)
        with open(args.chrom_colors) as fh:
            chrom_colors = {}
            for line in fh:
                chrom, color = line.strip().split("\t")
                if chrom not in _chroms:
                    continue
                chrom_colors[chrom] = color
    else:
        chrom_colors = DEF_CHR_COLORS

    if "all" not in chrom_colors:
        chrom_colors["all"] = "#808080"

    df_lengths = pl.read_csv(
        args.infile,
        has_header=False,
        separator="\t",
        columns=[0, 1, 2, 3],
        new_columns=DEF_COLS,
    ).with_columns(source=pl.lit("samples"))
    if plot_all:
        df_lengths = df_lengths.with_columns(chrom_name=pl.lit("all"))
    else:
        df_lengths = (
            df_lengths.with_columns(
                mtch_chrom=pl.col("chrom").str.extract_groups(rgx_name_groups),
            )
            .unnest("mtch_chrom")
            .drop("2")
        )

    added_palettes = OrderedDict()
    dfs_added_lengths = []
    if args.added_inputs:
        for elems in (lbl_path.split("=") for lbl_path in args.added_inputs):
            # Allow color to be optional.
            try:
                lbl, path, color = elems
            except ValueError:
                lbl, path = elems
                # Generate random color.
                color = np.random.rand(3)

            df = pl.read_csv(
                path,
                has_header=False,
                separator="\t",
                columns=[0, 1, 2, 3],
                new_columns=DEF_COLS,
            ).with_columns(source=pl.lit(lbl))

            if plot_all:
                df = df.with_columns(chrom_name=pl.lit("all"))
            else:
                df = df.with_columns(
                    chrom_name=pl.col("chrom").str.extract(f"^({rgx_chrom})$")
                )

            added_palettes[lbl] = color
            dfs_added_lengths.append(df)

    df_all_lengths: pl.DataFrame = pl.concat([df_lengths, *dfs_added_lengths])

    # Merge asat HOR array lengths
    if args.mode == "total":
        df_all_lengths = df_all_lengths.group_by(["chrom", "source"]).agg(
            pl.col("chrom_name").first(),
            pl.col("chrom_st").min(),
            pl.col("chrom_end").max(),
            pl.col("length").sum(),
        )
        ylabel = "Cumulative length of α-satellite HOR arrays (Mbp)"
    else:
        ylabel = "Length of α-satellite HOR arrays (Mbp)"

    df_all_lengths = df_all_lengths.with_columns(
        color_key=pl.when(pl.col("source") != "samples")
        .then(pl.col("source"))
        .otherwise(pl.col("chrom_name"))
    )
    palettes = chrom_colors | added_palettes
    # Add remaining chroms to plot if multi-chroms.
    uncovered_chroms = set(df_all_lengths["chrom_name"].unique()).difference(
        palettes.keys()
    )
    palettes = palettes | {chrom: "#FFFFFF" for chrom in uncovered_chroms}

    df_all_lengths_pd = df_all_lengths.to_pandas()
    sns.violinplot(
        x="chrom_name",
        y="length",
        hue="chrom_name",
        data=df_all_lengths_pd,
        palette=palettes,
        inner="quart",
        cut=0.75,
    )
    sns.swarmplot(
        x="chrom_name",
        y="length",
        data=df_all_lengths_pd,
        hue="color_key",
        linewidth=0.5,
        edgecolor="black",
        palette=palettes,
        size=4,
    )

    ax = plt.gca()
    legend_kwargs = dict(
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        title="Source",
        alignment="left",
        frameon=False,
    )
    try:
        # Sort legend elements
        legend_elem_order = {
            elem: i for i, elem in enumerate([*args.chroms, *added_palettes.keys()])
        }
        handles_labels = ax.get_legend_handles_labels()
        handles, labels = zip(
            *sorted(zip(*handles_labels), key=lambda x: legend_elem_order[x[1]])
        )

        # Place outside of figure.
        ax.legend(handles, labels, **legend_kwargs)
    except ValueError:
        ax.legend(
            handles=[
                Patch(color=color, label=chrom) for chrom, color in chrom_colors.items()
            ],
            **legend_kwargs,
        )

    # Hide spines
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    # x-axis
    ax.set_xlabel("Chromosome")
    # Remove chr from x-ticks
    xtick_labels = [lbl.get_text().replace("chr", "") for lbl in ax.get_xticklabels()]
    ax.set_xticklabels(xtick_labels)

    # Set units of y-axis
    ax.yaxis.minorticks_on()
    # Remove sci notation
    ax.yaxis.set_major_formatter("plain")
    new_xtick_labels = []
    _, ymax = ax.get_ylim()
    yticks, yticklabels = ax.get_yticks(), ax.get_yticklabels()
    for txt in yticklabels:
        _, y = txt.get_position()
        # Convert units and round.
        new_y_txt = str(round(y / 1_000_000, 3))
        txt.set_text(new_y_txt)
        new_xtick_labels.append(txt)

    # Add line and ytick for mean length
    mean_length = df_all_lengths["length"].mean()
    ax.axhline(mean_length, linestyle="dotted", color="black")
    mean_length_label = str(round(mean_length / 1_000_000, 1))
    ax.set_yticks([*yticks, mean_length], [*new_xtick_labels, mean_length_label])

    ax.set_ylabel(ylabel)
    ax.set_ylim(0, ymax)

    plt.gcf().set_size_inches(20, 8)
    plt.savefig(args.output, dpi=600, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
