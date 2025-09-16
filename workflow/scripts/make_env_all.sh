#!/bin/bash

set -euo pipefail

conda-merge $(find workflow/ -wholename "*/env*/*.y*ml") | \
    grep -v "name:" > env_all.yaml
