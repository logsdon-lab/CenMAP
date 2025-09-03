import os
import sys
import yaml
import copy
import json
import argparse
import polars as pl

from typing import Any
from collections import defaultdict

from cenplot import plot_one_cen, read_one_cen_tracks


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
            df = pl.read_csv(file, separator="\t", has_header=False)
            for prt, df_chrom in df.partition_by(["column_1"], as_dict=True).items():
                chrom = prt[0]
                fpath = os.path.join(outdir, f"{chrom}_{dtype}.bed")
                df_chrom.write_csv(fpath, separator="\t", include_header=False)
                bed_files[chrom][dtype] = fpath

    tracks = []
    for chrom, dtype_bedfiles in bed_files.items():
        for trk in track_format["tracks"]:
            dtype = trk["path"]
            bed_file = dtype_bedfiles.get(dtype)
            if not isinstance(bed_file, str):
                print(f"No data for {dtype}.", file=sys.stderr)
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
    spacer_track = {
        "position": "relative",
        "proportion": 0.01,
        "type": "spacer",
    }
    legend_track = {
        "position": "relative",
        "proportion": 0.25,
        "type": "legend",
        "options": {"index": args.ref_ax_idx, "legend_ncols": 10},
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
        tracks, settings = read_one_cen_tracks(fh, chrom=None)
        for track in tracks.tracks:
            if not isinstance(track.data, pl.DataFrame):
                continue
            # Update legend title.
            chrom = track.data["chrom"].first()
            if hasattr(track.options, "legend_title"):
                track.options.legend_title = chrom
        _ = plot_one_cen(tracks.tracks, outdir, bname, settings)

    if not args.keep_tempfiles:
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


if __name__ == "__main__":
    raise SystemExit(main())
