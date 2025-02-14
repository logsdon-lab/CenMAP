import os
import sys
import copy
import argparse
import multiprocessing

import numpy as np
import polars as pl
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from typing import Iterable, TextIO, BinaryIO
from concurrent.futures import ProcessPoolExecutor

from cenplot import (
    PlotSettings,
    plot_one_cen,
    read_one_cen_tracks,
    Track,
)


def create_shared_legend(
    plots: list[tuple[Figure, np.ndarray, str]],
    outdir: str,
    reference_ax_idx: int,
    dpi: int | None = None,
    transparent: bool = True,
):
    legend_fig, legend_axes, _ = copy.deepcopy(plots[reference_ax_idx])
    legend_ax: Axes = legend_axes[0, 0]

    all_legend_elems = {}
    for _, axes, _ in plots:
        ax: Axes = axes[0, 0]
        # Legend items.
        legend_handles, legend_labels = ax.get_legend_handles_labels()
        legend_elems = dict(zip(legend_labels, legend_handles))
        all_legend_elems.update(legend_elems)

    # Clear copied axis.
    legend_fig.suptitle(None)
    for ax in legend_fig.axes:
        ax.clear()
        # Remove all patches.
        for ptch in ax.patches:
            ptch.set_visible(False)
        # Remove all ticks and bottom spine.
        ax.set_yticks([], [])
        ax.set_xticks([], [])
        ax.axes.spines["bottom"].set_visible(False)

    legend = legend_ax.legend(
        legend_elems.values(),
        legend_elems.keys(),
        handlelength=1.0,
        handleheight=1.0,
        ncols=10,
        frameon=False,
        loc="center",
        alignment="center",
    )
    # Set patches edge color manually.
    # Turns out get_legend_handles_labels will get all rect patches and setting linewidth will cause all patches to be black.
    for ptch in legend.get_patches():
        ptch.set_linewidth(1.0)
        ptch.set_edgecolor("black")

    legend_outfile = os.path.join(outdir, "legend.png")
    legend_fig.tight_layout()
    legend_fig.savefig(legend_outfile, dpi=dpi, transparent=transparent)

    return (legend_fig, legend_axes, [legend_outfile])


def main():
    ap = argparse.ArgumentParser(
        description="Draw centromere tracks.",
    )
    ap.add_argument(
        "-t",
        "--input_tracks",
        required=True,
        type=argparse.FileType("rb"),
        help=(
            "TOML file with headerless BED files to plot. "
            "Specify under tracks the following fields: {name, position, type, proportion, path, or options}."
        ),
    )
    ap.add_argument(
        "-c",
        "--chroms",
        type=argparse.FileType("rt"),
        help="Names to plot in this order. Corresponds to 1st col in BED files.",
        required=True,
    )
    ap.add_argument(
        "-d",
        "--outdir",
        help="Output dir to plot multiple separate figures.",
        type=str,
        default=".",
    )
    ap.add_argument("--share_xlim", help="Share x-axis limits.", action="store_true")
    ap.add_argument(
        "--ref_ax_idx",
        help="Index of reference plot to add legend.",
        default=-1,
        type=int,
    )
    ap.add_argument("-p", "--processes", type=int, default=4, help="Processes to run.")
    args = ap.parse_args()

    input_tracks: BinaryIO = args.input_tracks
    chroms: TextIO = args.chroms
    outdir: str = args.outdir
    share_xlim: bool = args.share_xlim
    processes: int = args.processes
    all_chroms: Iterable[str] = [line.strip() for line in chroms.readlines()]

    inputs: list[tuple[Track, str, str, PlotSettings]] = []
    tracks_settings = []
    for chrom in all_chroms:
        try:
            tracks_settings.append(
                (chrom, *read_one_cen_tracks(input_tracks, chrom=chrom))
            )
        except Exception:
            pass

    dpi = None
    transparent = True
    xmin_all, xmax_all = sys.maxsize, 0
    if share_xlim:
        for *_, settings in tracks_settings:
            if settings.xlim:
                xmin, xmax = settings.xlim
                xmin_all = min(xmin_all, xmin)
                xmax_all = max(xmax_all, xmax)
            if settings.dpi:
                dpi = settings.dpi
            if settings.transparent:
                transparent = settings.transparent

    # Add position and legend track at end.
    for chrom, tracks_summary, plot_settings in tracks_settings:
        if share_xlim:
            plot_settings.xlim = (xmin_all, xmax_all)
        tracks = [
            Track(
                trk.title,
                trk.pos,
                trk.opt,
                trk.prop,
                trk.data.filter(pl.col("chrom") == chrom)
                if isinstance(trk.data, pl.DataFrame)
                else None,
                trk.options,
            )
            for trk in tracks_summary.tracks
        ]
        inputs.append(
            (
                tracks,
                outdir,
                chrom,
                plot_settings,
            )
        )

    os.makedirs(outdir, exist_ok=True)
    if processes == 1:
        plots = [plot_one_cen(*draw_arg) for draw_arg in inputs]
    else:
        with ProcessPoolExecutor(
            max_workers=processes, mp_context=multiprocessing.get_context("spawn")
        ) as pool:
            futures = [
                (draw_arg[2], pool.submit(plot_one_cen, *draw_arg))
                for draw_arg in inputs
            ]  # type: ignore[assignment]
            plots = []
            for chrom, future in futures:
                if future.exception():
                    print(
                        f"Failed to plot {chrom} ({future.exception()})",
                        file=sys.stderr,
                    )
                    continue
                plots.append(future.result())

    legend_plot = create_shared_legend(
        plots, outdir, args.ref_ax_idx, dpi=dpi, transparent=transparent
    )
    plots.append(legend_plot)

    merged_images = np.concatenate([plt.imread(files[0]) for _, _, files in plots])
    plt.imsave(f"{outdir}.png", merged_images)
    plt.imsave(f"{outdir}.pdf", merged_images)


if __name__ == "__main__":
    raise SystemExit(main())
