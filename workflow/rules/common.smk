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


# Include contants so can run individual smk files without main Snakefile.
include: "constants.smk"
