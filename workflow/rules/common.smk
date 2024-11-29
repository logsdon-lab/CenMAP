import os
import re
import sys
import json
import glob
from collections import defaultdict
from snakemake.settings.types import DeploymentMethod


# Include contants so can run individual smk files without main Snakefile.
include: "constants.smk"


# TODO: Migrate to YAML and use Snakemake built-in config validators.
try:
    with open(config["dna_brnn"]["full_alr_thr_file"]) as fh:
        DNA_BRNN_FULL_ALR_THRS = json.load(fh)
        DNA_BRNN_DEF_FULL_ALR_THR = DNA_BRNN_FULL_ALR_THRS.get(
            "default", DEF_DNA_BRNN_FULL_ALR_THR
        )
except (KeyError, FileNotFoundError):
    DNA_BRNN_FULL_ALR_THRS = {c: DEF_DNA_BRNN_FULL_ALR_THR for c in CHROMOSOMES}
except json.decoder.JSONDecodeError as json_err:
    raise Exception(
        f"Invalid JSON configuration for {config['dna_brnn']['full_alr_thr_file']}: {json_err}"
    )

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


def dnabrnn_alr_region_threshold(wc) -> int:
    return DNA_BRNN_FULL_ALR_THRS.get(str(wc.chr), DNA_BRNN_DEF_FULL_ALR_THR)


def censtats_status_edge_len(wc) -> int:
    return CENSTATS_STATUS_FULL_EDGE_LEN_THR.get(
        str(wc.chr), DEF_CENSTATS_STATUS_EDGE_LEN_THR
    )


def censtats_status_edge_perc_alr_thr(wc) -> int:
    return CENSTATS_STATUS_FULL_EDGE_PERC_ALR_THR.get(
        str(wc.chr), DEF_CENSTATS_STATUS_EDGE_PERC_ALR_THR
    )


def censtats_status_max_alr_len_thr(wc) -> int:
    return CENSTATS_STATUS_FULL_MAX_ALR_LEN_THR.get(
        str(wc.chr), DEF_CENSTATS_STATUS_MAX_ALR_LEN_THR
    )


def extract_fnames_and_chr(
    pattern: str, *, filter_chr: str | None = None
) -> tuple[list[str], list[str]]:
    fnames = glob_wildcards(pattern).fname
    filtered_fnames, chrs = [], []
    for fname in fnames:
        if mtch_chr_name := re.search(RGX_CHR, fname):
            chr_name = mtch_chr_name.group().strip("_")

            if not filter_chr:
                filtered_fnames.append(fname)
                chrs.append(chr_name)
                continue

            # Filter by chr.
            if chr_name != filter_chr:
                continue

            filtered_fnames.append(fname)
            chrs.append(chr_name)

    assert len(filtered_fnames) == len(
        chrs
    ), f"One or more fa files ({pattern}) does not contain a chromosome in its name."

    return filtered_fnames, chrs
