import os
import sys
import yaml
import copy
import json
import argparse
import polars as pl

from typing import Any
from collections import defaultdict

from cenplot import plot_tracks, read_tracks


def cleanup(cfg: str, bed_files: defaultdict[str, dict[str, str]]):
    try:
        os.remove(cfg)
    except FileNotFoundError:
        pass
    for chrom, dtype_bedfiles in bed_files.items():
        for dtype, file in dtype_bedfiles.items():
            try:
                os.remove(file)
            except FileNotFoundError:
                pass


def main():
    ap = argparse.ArgumentParser(
        description="Draw centromere tracks for multiple centromeres.",
    )
    ap.add_argument(
        "-i", "--infiles", required=True, type=str, help="Infile JSON files by type."
    )
    ap.add_argument(
        "-t",
        "--track_format",
        required=True,
        type=argparse.FileType("rb"),
        help=(
            "TOML file for centromere plot. "
            "Each track's path should have format string matching infile JSON str key."
        ),
    )
    ap.add_argument(
        "-o",
        "--output_prefix",
        help="Output prefix to plot multiple separate figures.",
        type=str,
        default="out",
    )
    ap.add_argument(
        "--ref_ax_idx",
        help="Index of reference plot to add legend.",
        default=0,
        type=int,
    )
    ap.add_argument(
        "--options",
        type=str,
        help="Options by dtype to replace options. ex. cfg['hor']['live_only']",
        default=None,
    )
    ap.add_argument(
        "--omit_if_any_empty",
        action="store_true",
        help="Omit track if any file is empty.",
    )
    ap.add_argument(
        "--keep_tempfiles",
        action="store_true",
        help="Keep files used to generate config.",
    )
    args = ap.parse_args()

    infiles: dict[str, list[str]] = json.loads(args.infiles)
    track_format: dict[str, Any] = yaml.safe_load(args.track_format)
    output_prefix: str = args.output_prefix
    outdir = os.path.dirname(output_prefix)
    outdir = os.path.abspath(outdir if outdir else ".")
    bname = os.path.basename(output_prefix)
    options: dict[str, dict[str, Any]] = (
        yaml.safe_load(args.options) if args.options else {}
    )

    # Read in all files by dtype.
    bed_files: defaultdict[str, dict[str, str]] = defaultdict(dict)

    for dtype, files in infiles.items():
        if isinstance(files, str):
            files = [files]
        elif isinstance(files, list):
            files = files
        else:
            raise ValueError(f"Unexpected type for {files}")
        for file in files:
            try:
                df = pl.read_csv(file, separator="\t", has_header=False)
            except pl.exceptions.NoDataError:
                continue
            for prt, df_chrom in df.partition_by(["column_1"], as_dict=True).items():
                chrom = prt[0]
                fpath = os.path.join(outdir, f"{chrom}_{dtype}.bed")
                df_chrom.write_csv(fpath, separator="\t", include_header=False)
                bed_files[chrom][dtype] = fpath

    dtypes = set(infiles.keys())
    spacer_track = {
        "position": "relative",
        "proportion": 0.01,
        "type": "spacer",
    }

    idx = 0
    tracks = []
    ref_indices = []
    missing_chroms = set()
    for chrom, dtype_bedfiles in bed_files.items():
        for i, trk in enumerate(track_format["tracks"]):
            dtype = trk["path"]
            bed_file = dtype_bedfiles.get(dtype)
            takes_space = trk.get("proportion")
            has_data = isinstance(bed_file, str)

            if takes_space:
                if i == args.ref_ax_idx:
                    ref_indices.append(idx)
                idx += 1

            if not has_data:
                print(f"No data for {chrom} {dtype}.", file=sys.stderr)
                # Add spacer if data not present.
                if takes_space:
                    spacer = copy.deepcopy(spacer_track)
                    spacer["proportion"] = trk["proportion"]
                    tracks.append(spacer)

                # Has no data and is an input file.
                if dtype in dtypes:
                    missing_chroms.add(chrom)
                continue

            new_trk = copy.deepcopy(trk)
            # Update config.
            added_options = options.get(trk["type"])
            if added_options:
                for k, v in added_options.items():
                    new_trk["options"][k] = v

            new_trk["path"] = bed_file
            tracks.append(new_trk)

    position_track = {
        "position": "relative",
        "proportion": 0.05,
        "type": "position",
        "options": {"hide_x": False},
    }
    legend_track = {
        "position": "relative",
        "proportion": 0.25,
        "type": "legend",
        "options": {"index": ref_indices, "legend_ncols": 10},
    }
    tracks.append(position_track)
    tracks.append(copy.deepcopy(spacer_track))
    tracks.append(legend_track)
    tracks.append(copy.deepcopy(spacer_track))

    plot_settings = track_format.get("settings")
    track_format["settings"]["dim"] = [
        20,
        (plot_settings["dim"][1] * len(bed_files.keys())) + 2,
    ]
    track_format["tracks"] = tracks
    cfg = os.path.join(f"{output_prefix}.yaml")
    with open(cfg, "wt") as fh:
        yaml.safe_dump(track_format, fh, allow_unicode=True)

    with open(cfg, "rb") as fh:
        try:
            tracks, settings = read_tracks(fh)
        except Exception as err:
            os.remove(cfg)
            raise RuntimeError(f"Failed to read tracks. {err}")

    final_tracks = []
    for track in tracks.tracks:
        if not track.data.is_empty():
            # Update legend title.
            chrom = track.data["chrom"].first()
            if chrom in missing_chroms and args.omit_if_any_empty:
                continue

            if hasattr(track.options, "legend_title"):
                track.options.legend_title = chrom

        final_tracks.append(track)

    if not final_tracks:
        os.remove(cfg)
        raise RuntimeError("No tracks to plot.")

    _ = plot_tracks(tracks=final_tracks, settings=settings, outdir=outdir, chrom=bname)

    if not args.keep_tempfiles:
        cleanup(cfg, bed_files)


if __name__ == "__main__":
    raise SystemExit(main())
