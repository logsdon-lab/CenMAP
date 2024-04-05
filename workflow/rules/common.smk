import os
import glob
import yaml
import pandas as pd
from collections import defaultdict


def get_hifi_read_wildcards() -> dict[str, list[str]]:
    """
    Get hifi reads by sample automatically from hifi_reads_dir.
    Expects {hifi_reads_dir}/{sample}/*.bam
    """
    reads_dir = config["nuc_freq"]["hifi_reads_dir"]
    # Avoid subdirs by constraining wildcards.
    # https://stackoverflow.com/a/60744040
    path_pattern = os.path.join(reads_dir, "{sm,[^/]+}", "{flowcell_id,[^/]+}.bam")
    reads_run_mdata_id = glob_wildcards(path_pattern)

    samples = defaultdict(list)
    for sm, flowcell_id in zip(reads_run_mdata_id.sm, reads_run_mdata_id.flowcell_id):
        samples[sm].append(flowcell_id)
    return samples


def build_awk_cen_region_length_thr(chr_name: str) -> str:
    region_thresholds: list[list[int | None, int | None]] = (
        config["dna_brnn"]
        .get("repeat_len_thr", {})
        .get(
            chr_name, [config["dna_brnn"].get("default_repeat_len_thr", [1_000, None])]
        )
    )
    stmts = []
    for thr_min, thr_max in region_thresholds:
        if thr_min and thr_max:
            stmts.append(f"$5>{thr_min} && $5<{thr_max}")
        elif thr_min:
            stmts.append(f"$5>{thr_min}")
        else:
            stmts.append(f"$5<{thr_max}")

    return f"({' || '.join(stmts)})"


# Include contants so can run individual smk files without main Snakefile.
include: "constants.smk"
