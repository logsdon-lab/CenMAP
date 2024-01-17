#!/bin/bash

set -euo pipefail 

# Either globus login or set register app and get tokens.
# https://docs.globus.org/api/auth/developer-guide/#developing-apps

if [[ $# -ne 2 ]]
then
    echo "usage: $0 asm_output_dir raw_data_output_dir" 1>&2
    exit 1
fi

ASM_DIR=$1
RAW_DATA_DIR=$2

if [ ! -d $ASM_DIR ] || [ ! -d $RAW_DATA_DIR ]; then
    echo "One or more directories are invalid." 1>&2
    exit 1
fi

# EMBL-EBI Public Data ID
ENDPT="47772002-3e5b-4fd3-b97c-18cee38d6df2"
LOCAL_ENDPT=$(globus endpoint local-id)

HGSVC3_WORKING_PATH="/1000g/ftp/data_collections/HGSVC3/working/"
ASM_BATCH_1_PATH="$HGSVC3_WORKING_PATH/20230818_verkko_batch1/assemblies"
ASM_BATCH_2_PATH="$HGSVC3_WORKING_PATH/20230927_verkko_batch2/assemblies"

# transfer assemblies to local endpt
# only hap1, hap2, and unassigned
# https://docs.globus.org/cli/reference/transfer/#examples
globus transfer \
    -r $ENDPT:$ASM_BATCH_1_PATH \
    $LOCAL_ENDPT:$ASM_DIR \
    --exclude "*contaminants*" \
    --exclude "*mito*" \
    --exclude "*rdna*"

globus transfer \
    -r $ENDPT:$ASM_BATCH_2_PATH \
    $LOCAL_ENDPT:$ASM_DIR \
    --exclude "*contaminants*" \
    --exclude "*mito*" \
    --exclude "*rdna*"

# Download raw ONT/HIFI.
globus ls $ENDPT:$HGSVC3_WORKING_PATH | \
    egrep -E "HiFi|HIFI|ONT" \
    xargs -P 8 -I [] \
    globus transfer -r "$ENDPT:$HGSVC3_WORKING_PATH/[]" $LOCAL_ENDPT:$RAW_DATA_DIR
