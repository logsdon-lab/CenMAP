import os
import re
import json
import glob
from collections import defaultdict


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
        CENSTATS_STATUS_DEF_EDGE_LEN_THR = CENSTATS_STATUS_FULL_EDGE_LEN_THR.get(
            "default", DEF_CENSTATS_STATUS_EDGE_LEN_THR
        )
        # Edge percent alpha-satellite repeat threshold to be considered partial.
        CENSTATS_STATUS_FULL_EDGE_PERC_ALR_THR = CENSTATS_STATUS_CFG[
            "edge_perc_alr_thr"
        ]
        CENSTATS_STATUS_DEF_EDGE_PERC_ALR_THR = (
            CENSTATS_STATUS_FULL_EDGE_PERC_ALR_THR.get(
                "default", DEF_CENSTATS_STATUS_EDGE_PERC_ALR_THR
            )
        )
        # Max required alpha-satellite repeat length threshold
        # Smaller thresholds are edge-case for chrs whose repeats are small and broken up.
        CENSTATS_STATUS_FULL_MAX_ALR_LEN_THR = CENSTATS_STATUS_CFG["max_alr_len_thr"]
        CENSTATS_STATUS_DEF_MAX_ALR_LEN_THR = CENSTATS_STATUS_FULL_MAX_ALR_LEN_THR.get(
            "default", DEF_CENSTATS_STATUS_MAX_ALR_LEN_THR
        )
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
        str(wc.chr), CENSTATS_STATUS_DEF_EDGE_LEN_THR
    )


def censtats_status_edge_perc_alr_thr(wc) -> int:
    return CENSTATS_STATUS_FULL_EDGE_PERC_ALR_THR.get(
        str(wc.chr), CENSTATS_STATUS_DEF_EDGE_PERC_ALR_THR
    )


def censtats_status_max_alr_len_thr(wc) -> int:
    return CENSTATS_STATUS_FULL_MAX_ALR_LEN_THR.get(
        str(wc.chr), CENSTATS_STATUS_DEF_MAX_ALR_LEN_THR
    )


def extract_fa_fnames_and_chr(input_dir: str) -> tuple[list[str], list[str]]:
    fnames = glob_wildcards(os.path.join(input_dir, "{fname}.fa")).fname
    chrs = []
    for fname in fnames:
        chr_name = re.search(RGX_CHR, fname)
        if chr_name:
            chrs.append(chr_name.group(0))

    assert len(fnames) == len(
        chrs
    ), f"One or more fa files in {input_dir} does not contain a chromosome in its name."

    return fnames, chrs
