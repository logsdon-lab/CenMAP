import sys
import argparse
import polars as pl
import cenplot as cplt
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

BED9_COLS = (
    "chrom",
    "chrom_st",
    "chrom_end",
    "name",
    "score",
    "strand",
    "chrom_tst",
    "chrom_tend",
    "item_rgb",
)


def format_xaxis_ticklabels(ax: Axes):
    # Remove scientific notation.
    ax.xaxis.set_major_formatter("plain")
    new_xtick_labels = []
    xmin, xmax = ax.get_xlim()
    xticks, xticklabels = list(ax.get_xticks()), ax.get_xticklabels()
    for txt in xticklabels:
        x, _ = txt.get_position()
        # Convert units and round.
        new_x_txt = round(x / 1_000_000, 1)
        txt.set_text(new_x_txt)
        new_xtick_labels.append(txt)

    # Add last position.
    xticks.append(xmin)
    xticks.append(xmax)
    new_xtick_labels.append("0.0")
    new_xtick_labels.append(str(round(xmax / 1_000_000, 1)))

    ax.set_xticks(xticks, new_xtick_labels, fontsize="medium")
    ax.set_xlabel("Position (Mbp)", fontsize="medium")
    ax.set_xlim(xmin, xmax)


def main():
    ap = argparse.ArgumentParser(description="Generate ideogram.")
    ap.add_argument(
        "-i",
        "--bed_cen",
        type=str,
        required=True,
        help="Complete and accurate centromere BED9 file.",
    )
    ap.add_argument(
        "-g",
        "--genome_sizes",
        type=str,
        required=True,
        help="Genome assembly contig lengths. A fasta index.",
    )
    ap.add_argument("-o", "--output", default="centromeres.pdf", help="Output plot.")
    ap.add_argument(
        "--color_cen", type=str, default="#ff0000", help="Color of centromere."
    )
    ap.add_argument("--color_ctg", type=str, default="#d3d3d3", help="Color of contig.")
    ap.add_argument("-t", "--title", type=str, default=None, help="Title of plot.")
    ap.add_argument(
        "-ht",
        "--height",
        type=int,
        default=20,
        help="Height of plot",
    )
    ap.add_argument(
        "-w",
        "--width",
        type=int,
        default=12,
        help="Width of plot",
    )
    ap.add_argument("--legend_prop", type=float, default=0.1, help="Legend proportion.")
    ap.add_argument(
        "--use_renamed_reoriented",
        action="store_true",
        help="Reorient relative to reference.",
    )
    args = ap.parse_args()

    try:
        df_cens = pl.read_csv(
            args.bed_cen, new_columns=BED9_COLS, separator="\t", has_header=False
        )
    except pl.exceptions.NoDataError:
        print("No complete cens.", file=sys.stderr)
        fig, axes = plt.subplots()
        fig.savefig(args.output, bbox_inches="tight", dpi=600, transparent=True)
        return 0

    genome_sizes = dict(
        pl.read_csv(
            args.genome_sizes,
            new_columns=["contig", "contig_length"],
            separator="\t",
            columns=[0, 1],
            has_header=False,
        ).iter_rows()
    )

    tracks = []
    max_contig_length = 0
    for grp, df_grp in df_cens.group_by(["chrom"]):
        contig = grp[0]
        contig_length = genome_sizes.get(contig)
        if not contig_length:
            print(f"Missing {contig} in genome sizes. Skipping.", file=sys.stderr)
            continue

        if args.use_renamed_reoriented:
            added_cmds = {
                "chrom_st": pl.col("chrom_tst"),
                "chrom_end": pl.col("chrom_tend"),
            }
            contig = df_grp["name"].first()
        else:
            added_cmds = {}

        # Set color to red.
        # item_rgb doesn't get set unless map_value_colors is called.
        df_grp = df_grp.with_columns(color=pl.lit(args.color_cen), **added_cmds)

        # Make row for whole contig.
        row_ctg = df_grp.row(0, named=True)
        row_ctg["chrom_st"] = 1
        row_ctg["chrom_tst"] = 1
        row_ctg["chrom_end"] = contig_length
        row_ctg["chrom_tend"] = contig_length
        row_ctg["color"] = args.color_ctg

        # Create df with cen row and ctg row.
        df_all = pl.concat(
            [
                pl.DataFrame([row_ctg]),
                df_grp,
            ]
        )
        track = cplt.Track(
            title=None,
            pos=cplt.TrackPosition.Relative,
            opt=cplt.TrackType.Label,
            prop=0.5,
            data=df_all,
            options=cplt.LabelTrackSettings(
                legend=True,
                legend_title_only=True,
                legend_title_fontsize="small",
                legend_title=contig,
                bg_border=True,
                edgecolor="black",
                shape="rect",
            ),
        )
        tracks.append((track, contig_length))
        max_contig_length = max(max_contig_length, contig_length)

    # Sort by contig length with largest at top.
    tracks.sort(key=lambda x: x[1], reverse=True)
    tracks = [trk[0] for trk in tracks]

    # Add spacer and position
    tracks.insert(
        0,
        cplt.Track(
            title=None,
            pos=cplt.TrackPosition.Relative,
            opt=cplt.TrackType.Spacer,
            prop=0.5,
            data=pl.DataFrame(),
            options=cplt.SpacerTrackSettings(),
        ),
    )
    tracks.append(
        cplt.Track(
            title=None,
            pos=cplt.TrackPosition.Relative,
            opt=cplt.TrackType.Position,
            data=pl.DataFrame(),
            prop=0.2,
            options=cplt.PositionTrackSettings(legend=False),
        )
    )
    idx_position = len(tracks)

    fig, axes, _ = cplt.plot_tracks(
        tracks=tracks,
        settings=cplt.PlotSettings(
            title=args.title,
            format="pdf",
            dim=(args.width, args.height),
            layout="tight",
            legend_pos=cplt.LegendPosition.Left,
            legend_prop=args.legend_prop,
            xlim=(1, max_contig_length),
            axis_h_pad=0.5,
        ),
    )
    format_xaxis_ticklabels(axes[idx_position - 1, 1])
    fig.savefig(args.output, bbox_inches="tight", dpi=600, transparent=True)


if __name__ == "__main__":
    raise SystemExit(main())
