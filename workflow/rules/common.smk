import re
import sys
import json
import glob
import tomllib
import yaml

from os.path import join, dirname, basename, splitext
from collections import defaultdict, Counter
from typing import Any
from snakemake.settings.types import DeploymentMethod
from snakemake.utils import validate


# Validate and fill defaults.
validate(config, "../../config/config.schema.yaml")


# Include contants so can run individual smk files without main Snakefile.
include: "constants.smk"
