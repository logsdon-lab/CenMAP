import os
import re
import json
import glob
from collections import defaultdict


with open(config["dna_brnn"]["full_alr_thr_file"]) as fh:
    DNA_BRNN_FULL_ALR_THRS = json.load(fh)
    DNA_BRNN_DEF_FULL_ALR_THR = DNA_BRNN_FULL_ALR_THRS.get("default", 1_000_000)


def get_hifi_read_wildcards() -> dict[str, list[str]]:
    """
    Get hifi reads by sample automatically from hifi_reads_dir.
    Expects {hifi_reads_dir}/{sample}/*.{ext}
    """
    # Avoid subdirs by constraining wildcards.
    # https://stackoverflow.com/a/60744040
    escaped_ext = re.escape("." + config["nuc_freq"].get("reads_ext", "bam"))
    path_pattern = re.compile(r"([^/]+)(" + escaped_ext + ")")
    samples = defaultdict(list)
    for root, read_dirs, _ in os.walk(config["nuc_freq"]["hifi_reads_dir"]):
        for read_dir in read_dirs:
            read_dir_path = os.path.join(root, read_dir)
            for file in os.listdir(read_dir_path):
                try:
                    flowcell_id, _ = re.search(path_pattern, file).groups()
                except (ValueError, AttributeError):
                    continue

                samples[read_dir].append(flowcell_id)

    return samples


def alr_region_threshold(wc) -> int:
    return DNA_BRNN_FULL_ALR_THRS.get(str(wc.chr), DNA_BRNN_DEF_FULL_ALR_THR)


# Include contants so can run individual smk files without main Snakefile.
include: "constants.smk"
