import os
import sys
import yaml
import tomllib
import argparse
from io import BytesIO
from typing import Any, TextIO


def format_toml_path(
    input_plot_layout: BytesIO,
    output_plot_layout: TextIO,
    options: dict[str, Any],
    path_kwargs: Any,
):
    """
    Format the plot layout filling any kwargs in the path.
    Also checks if path is empty.
    """
    settings = tomllib.load(input_plot_layout)
    new_settings = {"settings": settings["settings"], "tracks": []}

    number_skipped = 0
    for trk in settings["tracks"]:
        path = trk.get("path")
        new_trk = trk.copy()

        try:
            for opt, _ in trk.get("options", {}).items():
                new_value = options.get(opt)
                if new_value:
                    new_trk["options"][opt] = new_value
        except KeyError:
            print(
                f"Invalid format key in path {path} from cenplot track file, {input_plot_layout.name}.",
                file=sys.stderr,
            )
        # Some tracks don't have paths.
        if not path:
            # If any track is a legend track, need to track index and update when empty (-number skipped.)
            if trk["type"] == "legend":
                original_idx = trk.get("options", {}).get("index")
                if original_idx:
                    new_trk["options"]["index"] = max(0, original_idx - number_skipped)
            else:
                new_trk = trk

            new_settings["tracks"].append(new_trk)
            continue
        try:
            new_path = path.format(**path_kwargs)
        except KeyError:
            print(
                f"Invalid format key in path {path} from cenplot track file, {input_plot_layout.name}.",
                file=sys.stderr,
            )
            continue

        # Skip if empty.
        if os.stat(new_path).st_size == 0:
            number_skipped += 1
            continue

        new_trk = trk.copy()
        # Pass params from snakemake
        new_trk["path"] = new_path
        new_settings["tracks"].append(new_trk)

    output_plot_layout.write(yaml.dump(new_settings))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--infile", type=argparse.FileType("rb"), required=True)
    ap.add_argument("-o", "--outfile", type=argparse.FileType("wt"), default=sys.stdout)
    ap.add_argument(
        "-p",
        "--path_kwargs",
        nargs="*",
        metavar="{key}={value}",
        type=lambda x: x.split("=", 1),
    )
    ap.add_argument(
        "--options",
        nargs="*",
        metavar="{key}={value}",
        type=lambda x: x.split("=", 1),
    )
    args = ap.parse_args()

    format_toml_path(
        input_plot_layout=args.infile,
        output_plot_layout=args.outfile,
        path_kwargs=dict(args.path_kwargs) if args.path_kwargs else {},
        options=dict(args.options) if args.options else {},
    )


if __name__ == "__main__":
    raise SystemExit(main())
