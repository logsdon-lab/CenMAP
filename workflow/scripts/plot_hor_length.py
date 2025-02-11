import argparse
import numpy as np
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt

from collections import OrderedDict


DEF_COLS = ["chrom", "chrom_st", "chrom_end", "length"]
CHR_ORDER = [f"chr{i}" for i in (*range(1, 23), "X", "Y")]
CHR_COLORS = dict(
    zip(
        CHR_ORDER,
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
        help="Additional centromere HOR array lengths.",
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
    ap.add_argument("-o", "--output", help="Output plot file.", required=True, type=str)

    args = ap.parse_args()

    df_lengths = pl.read_csv(
        args.infile,
        has_header=False,
        separator="\t",
        columns=[0, 1, 2, 3],
        new_columns=DEF_COLS,
    ).with_columns(source=pl.lit("samples"))

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
            added_palettes[lbl] = color
            dfs_added_lengths.append(df)

    df_all_lengths: pl.DataFrame = pl.concat([df_lengths, *dfs_added_lengths])
    df_all_lengths = df_all_lengths.with_columns(
        chrom_name=pl.col("chrom")
        .str.extract(r"(chr\d+|chrX|chrY)")
        .cast(pl.Enum(CHR_ORDER))
    )

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

    palettes = CHR_COLORS | added_palettes
    sns.violinplot(
        x="chrom_name",
        y="length",
        hue="chrom_name",
        data=df_all_lengths,
        palette=palettes,
        inner="quart",
        cut=0.75,
    )
    sns.swarmplot(
        x="chrom_name",
        y="length",
        data=df_all_lengths,
        hue="color_key",
        linewidth=0.5,
        edgecolor="black",
        palette=palettes,
        size=4,
    )

    ax = plt.gca()

    # Sort legend elements
    legend_elem_order = {
        elem: i for i, elem in enumerate([*CHR_ORDER, *added_palettes.keys()])
    }
    handles_labels = ax.get_legend_handles_labels()
    handles, labels = zip(
        *sorted(zip(*handles_labels), key=lambda x: legend_elem_order[x[1]])
    )

    # Place outside of figure.
    ax.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        title="Source",
        alignment="left",
        frameon=False,
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
