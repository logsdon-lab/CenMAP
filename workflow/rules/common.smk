import re
import sys
import json
import glob
import tomllib
import yaml

from os.path import join
from collections import defaultdict, Counter
from typing import Any
from snakemake.settings.types import DeploymentMethod
from snakemake.utils import validate


# Validate and fill defaults.
validate(config, "../../config/config.schema.yaml")


# Include contants so can run individual smk files without main Snakefile.
include: "constants.smk"


# TODO: Migrate to YAML and use Snakemake built-in config validators.
try:
    with open(config["repeatmasker"]["config_censtats_status"]) as fh:
        CENSTATS_STATUS_CFG = json.load(fh)
        # Edge length to evaluate
        CENSTATS_STATUS_FULL_EDGE_LEN_THR = CENSTATS_STATUS_CFG["edge_len"]
        # Edge percent alpha-satellite repeat threshold to be considered partial.
        CENSTATS_STATUS_FULL_EDGE_PERC_ALR_THR = CENSTATS_STATUS_CFG[
            "edge_perc_alr_thr"
        ]
        # Max required alpha-satellite repeat length threshold
        # Smaller thresholds are edge-case for chrs whose repeats are small and broken up.
        CENSTATS_STATUS_FULL_MAX_ALR_LEN_THR = CENSTATS_STATUS_CFG["max_alr_len_thr"]
except (KeyError, FileNotFoundError):
    CENSTATS_STATUS_FULL_EDGE_LEN_THR = {
        c: DEF_CENSTATS_STATUS_EDGE_LEN_THR for c in CHROMOSOMES
    }
    CENSTATS_STATUS_FULL_EDGE_PERC_ALR_THR = {
        c: DEF_CENSTATS_STATUS_EDGE_PERC_ALR_THR for c in CHROMOSOMES
    }
    CENSTATS_STATUS_FULL_MAX_ALR_LEN_THR = {
        c: DEF_CENSTATS_STATUS_MAX_ALR_LEN_THR for c in CHROMOSOMES
    }
except json.decoder.JSONDecodeError as json_err:
    raise Exception(
        f"Invalid JSON configuration for {config['dna_brnn']['config_censtats_status_status']}: {json_err}"
    )


def censtats_status_edge_len(chrom: str) -> int:
    return CENSTATS_STATUS_FULL_EDGE_LEN_THR.get(
        chrom, DEF_CENSTATS_STATUS_EDGE_LEN_THR
    )


def censtats_status_edge_perc_alr_thr(chrom: str) -> int:
    return CENSTATS_STATUS_FULL_EDGE_PERC_ALR_THR.get(
        chrom, DEF_CENSTATS_STATUS_EDGE_PERC_ALR_THR
    )


def censtats_status_max_alr_len_thr(chrom: str) -> int:
    return CENSTATS_STATUS_FULL_MAX_ALR_LEN_THR.get(
        chrom, DEF_CENSTATS_STATUS_MAX_ALR_LEN_THR
    )


def get_chrom_name(name: str) -> str | None:
    if mtch_chr_name := re.search(RGX_CHR, name):
        return mtch_chr_name.group()


def format_toml_path(
    input_plot_layout: str, output_plot_layout: str, indir: str, **kwargs: Any
):
    """
    Format the plot layout filling in {indir} and any kwargs in the path.
    Also checks if path is empty.
    """
    # infile indicates where all other bedfiles are.
    with (
        open(input_plot_layout, "rb") as fh,
        open(output_plot_layout, "wt") as out_fh,
    ):
        settings = tomllib.load(fh)
        new_settings = {"settings": settings["settings"], "tracks": []}

        for trk in settings["tracks"]:
            path = trk.get("path")
            # Some tracks don't have paths.
            if not path:
                new_settings["tracks"].append(trk)
                continue
            try:
                new_path = path.format(indir=indir, **kwargs)
            except KeyError:
                print(
                    f"Invalid format key in path {path} from cenplot track file, {input.plot_layout}.",
                    file=sys.stderr,
                )
                continue

            # Skip if empty.
            if os.stat(new_path).st_size == 0:
                continue

            new_trk = trk.copy()
            # Pass params from snakemake
            new_trk["path"] = new_path
            new_settings["tracks"].append(new_trk)

        out_fh.write(yaml.dump(new_settings))
