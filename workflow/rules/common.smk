import os
import re
import json
import glob
from collections import defaultdict


with open(config["dna_brnn"]["full_alr_thr_file"]) as fh:
    DNA_BRNN_FULL_ALR_THRS = json.load(fh)
    DNA_BRNN_DEF_FULL_ALR_THR = DNA_BRNN_FULL_ALR_THRS.get("default", 1_000_000)


def alr_region_threshold(wc) -> int:
    return DNA_BRNN_FULL_ALR_THRS.get(str(wc.chr), DNA_BRNN_DEF_FULL_ALR_THR)


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
    assert len(chrs) > 0, f"No fa files found in {input_dir}."

    return fnames, chrs


# Include contants so can run individual smk files without main Snakefile.
include: "constants.smk"
