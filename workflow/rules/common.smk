import os
import glob
from collections import defaultdict


def get_hifi_read_wildcards() -> dict[str, list[str]]:
    """
    Get hifi reads by sample automatically from hifi_reads_dir.
    Expects {hifi_reads_dir}/{sample}/*.(bam|fq|fastq)(.gz)*
    """
    # Avoid subdirs by constraining wildcards.
    # https://stackoverflow.com/a/60744040
    path_pattern = re.compile(r"([^/]+)(\.bam|\.fq.gz|\.fastq\.gz|\.fq|\.fastq)")

    samples = defaultdict(list)
    for root, read_dirs, _ in os.walk(config["nuc_freq"]["hifi_reads_dir"]):
        for read_dir in read_dirs:
            read_dir_path = os.path.join(root, read_dir)
            for file in os.listdir(read_dir_path):
                try:
                    flowcell_id, extension = re.search(path_pattern, file).groups()
                except ValueError:
                    continue
                
                samples[read_dir].append(flowcell_id + extension)

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
