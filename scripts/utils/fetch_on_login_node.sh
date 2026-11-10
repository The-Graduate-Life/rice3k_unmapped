#!/bin/bash
# Fallback when compute nodes have no internet: fetch unmapped reads on a login/transfer node.
# Light on CPU (~7 s and ~190 MB per accession). Keep parallelism modest to be a good citizen.
#   bash scripts/utils/fetch_on_login_node.sh metadata/accessions_pilot100.txt 4
set -euo pipefail
LIST=${1:?list}; J=${2:-4}
cd "$(dirname "$0")/../.."; source config.sh; activate_env
xargs -a "$LIST" -P "$J" -I{} bash scripts/01_process_accession.sh {} 1 --fetch-only
du -sh "$WORK/fastq"
