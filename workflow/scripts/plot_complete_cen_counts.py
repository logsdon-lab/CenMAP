import argparse
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt

DEF_COLS = ["sample", "cnt", "perc"]
NUM_CHRS = 46


def main():
    ap = argparse.ArgumentParser(
        description="Plot centromere correct HOR array counts."
    )
    ap.add_argument(
        "-i",
        "--infile",
        help="Input correct centromere HOR array percentage and counts.",
        type=argparse.FileType("rb"),
    )
    ap.add_argument("-o", "--output", help="Output plot file.", type=str, required=True)

    args = ap.parse_args()

    df_cnts = pl.read_csv(
        args.infile, separator="\t", has_header=False, new_columns=DEF_COLS
    )
    df_plot = df_cnts.filter(pl.col("sample") != "all").sort(by="cnt", descending=True)
    sns.barplot(data=df_plot, x="sample", y="perc")
    (
        fig,
        ax,
    ) = plt.gcf(), plt.gca()
    ax.set_ylim(0, 100)

    # Mean
    mean_perc = round(df_plot["perc"].mean())
    ax.axhline(mean_perc, linestyle="dotted", color="black")
    yticks, ytick_labels = ax.get_yticks(), ax.get_yticklabels()
    # Add mean ytick.
    ax.set_yticks([*yticks, mean_perc], [*ytick_labels, str(mean_perc)])
    ax.set_ylabel(r"% of centromeres completely and accurately assembled")

    # Add secondary axis for percent correctly annotated.
    new_yticks = [*range(0, 50, 10), NUM_CHRS]
    ax_2 = ax.secondary_yaxis(location="right")
    ax_2.set_yticks(
        [round((ytick_lbl / NUM_CHRS) * 100) for ytick_lbl in new_yticks],
        [str(ytick_lbl) for ytick_lbl in new_yticks],
    )
    ax_2.set_ylabel("# of centromeres completely and accurately assembled")

    # Rotate x-axis labels for visibility
    ax.set_xlabel("Sample")
    # https://github.com/matplotlib/matplotlib/issues/13774#issuecomment-478250353
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # Hide spines
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    fig.set_size_inches(20, 8)
    plt.savefig(args.output, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
