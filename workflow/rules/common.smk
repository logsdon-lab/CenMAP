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
