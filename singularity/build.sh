#!/bin/bash
set -euo pipefail

# Builds the biopandas container for the
# prepare_reference and process_counts, pipeline steps

SCRIPT_DIR="$(dirname -- $(readlink -f $0))"
DEF_FILE="${SCRIPT_DIR}/biopandas.def"
IMG_FILE="${SCRIPT_DIR}/biopandas.sif"

unset SINGULARITY_BINDPATH
sudo -E singularity build --force $IMG_FILE $DEF_FILE
