import os
import sys
import yaml
import json
import argparse
from io import BytesIO
from typing import Any


def format_yaml_path(
    infiles: dict[str, str],
    input_plot_layout: BytesIO,
    options: dict[str, Any],
):
    """
    Format the plot layout filling any kwargs in the path.
    Also checks if path is empty.
    """
    settings: dict[str, Any] = yaml.safe_load(input_plot_layout)
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
            new_path = os.path.abspath(infiles[path])
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

    print(yaml.dump(new_settings))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i", "--infiles", type=str, required=True, help="Infiles as JSON map."
    )
    ap.add_argument("-t", "--track_format", type=argparse.FileType("rb"), required=True)
    ap.add_argument(
        "--options",
        type=str,
        help="Options by dtype to replace options. ex. cfg['hor']['live_only']",
        default=None,
    )
    args = ap.parse_args()

    infiles: dict[str, str] = json.loads(args.infiles)
    options: dict[str, dict[str, Any]] = (
        json.loads(args.options) if args.options else {}
    )
    format_yaml_path(
        infiles=infiles,
        input_plot_layout=args.track_format,
        options=options,
    )


if __name__ == "__main__":
    raise SystemExit(main())
