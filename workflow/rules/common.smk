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
    fname_elems = []
    for fasta in fastas:
        bname, _ = splitext(basename(fasta))
        fname, coord = bname.rsplit(":", 1)
        if filter_chrom:
            sm, chrom, ctg = RGX_SM_CHR_CTG.search(fname).groups()
        else:
            sm, ctg = RGX_SM_CTG.search(fname).groups()
            chrom = None
        if chrom:
            elems = (sm, chrom, ctg, coord)
        else:
            elems = (sm, ctg, coord)
        fname_elems.append(elems)

    fnames = []
    sort_key = lambda x: (x[1], x[2], x[3]) if filter_chrom else (x[1], x[2])
    # Sort by coords so if multiple chr, chr position in name (chr3-chr21) matches.
    sorted_wcs = sorted(
        fname_elems,
        key=sort_key,
        reverse=True,
    )

    # Store index of chrom per contig.
    # In cases of dicentric contigs (chr3-chr21). Don't want to include twice.
    ctg_counter = Counter()
    for elems in sorted_wcs:
        if filter_chrom:
            sm, chrom, ctg, coord = elems
            chrom_names: list[str] = chrom.replace("rc-", "").split("-")
        else:
            sm, ctg, coord = elems
            chrom = None
            chrom_names = []

        if filter_chrom and not filter_chrom in chrom_names:
            continue

        ctg_id = (sm, chrom, ctg, coord)
        idx = ctg_counter[ctg_id]
        ctg_counter[ctg_id] += 1
        if filter_chrom and filter_chrom != chrom_names[idx]:
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
