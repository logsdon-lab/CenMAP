import re
import sys
import json
import glob
import tomllib
import yaml

from os.path import join, dirname, basename, splitext
from collections import defaultdict, Counter
from typing import Any, NamedTuple
from snakemake.settings.types import DeploymentMethod
from snakemake.utils import validate
from snakemake.io import Wildcards

# Validate and fill defaults.
validate(config, "../../config/config.schema.yaml")


def chrom_coord_to_tsv(chrom_coord: str) -> str:
    chrom, coords = chrom_coord.rsplit(":", 1)
    st, end = coords.split("-")
    return f"{chrom}\\t{st}\\t{end}\\n"


class SplitFilename(NamedTuple):
    sm: str
    chrom: str | None
    ctg: str
    coord: str

    def __str__(self):
        if self.chrom:
            return f"{self.sm}_{self.chrom}_{self.ctg}:{self.coord}"
        else:
            return f"{self.sm}_{self.ctg}:{self.coord}"


def get_valid_fnames(
    fastas: list[str], *, filter_chrom: str | None
) -> list[SplitFilename]:
    """
    Extract wildcards and filter fasta list for chrom.

    # Args
    * fastas - list of fasta files. must have file extension. ex. {fname}.fa
    * filter_chrom - chrom to filter for.

    # Returns
    * list of SplitFilenames
    """
    fnames = []
    for fasta in fastas:
        bname, _ = splitext(basename(fasta))
        fname, coord = bname.rsplit(":", 1)
        if filter_chrom:
            sm, chrom, ctg = RGX_SM_CHR_CTG.search(fname).groups()
            chrom_names: list[str] = chrom.replace("rc-", "").split("-")
        else:
            sm, ctg = RGX_SM_CTG.search(fname).groups()
            chrom = None
            chrom_names = []

        if filter_chrom and not filter_chrom in chrom_names:
            continue

        fnames.append(SplitFilename(sm=sm, chrom=chrom, ctg=ctg, coord=coord))

    return fnames


def sort_chrom_order(chromosomes: list[str]) -> list[str]:
    if not chromosomes:
        return []

    chroms_sex = []
    for chrom in ("chrX", "chrY"):
        try:
            chromosomes.remove(chrom)
            chroms_sex.append(chrom)
        except ValueError:
            pass

    chroms_numbers = sorted(int(chrom.replace("chr", "")) for chrom in chromosomes)
    chroms_final = [f"chr{chrom}" for chrom in chroms_numbers]
    chroms_final.extend(chroms_sex)
    return chroms_final


# Include contants so can run individual smk files without main Snakefile.
include: "constants.smk"


def cmd_filter_fa_chrom(*added_cmds) -> str:
    if not CHROMOSOMES:
        return " | ".join(added_cmds)

    # Sort in reverse order so chr12 matched first and not chr1.
    sorted_chroms = sort_chrom_order(CHROMOSOMES.copy())
    sorted_chroms.reverse()

    rgx_chrom = "|".join(sorted_chroms)
    rgx_chrom = f"[_-]{rgx_chrom}[_-]"
    # "[_-](chr12|chr11|chr1)[_-]"
    cmd = f"seqkit grep -r -p '{rgx_chrom}'"
    for acmd in added_cmds:
        cmd += f" | {acmd}"
    return cmd
