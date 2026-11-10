#!/bin/bash
# Scheduler-independent job body; submitted by run_pipeline.sh (PBS, SLURM or local).
# Receives PROJ, WORK, LIST, NCPU (and TASK where noted) from the environment.
set -uo pipefail
cd "$PROJ" && source config.sh && activate_env || exit 1
mkdir -p "$PROJ/logs"
# Centred kinship matrix from the LD-pruned SNP set (GEMMA), computed once.
set -e
cd "$WORK/gwas"
[ -s output/kin.cXX.txt ] && { echo "kinship exists"; exit 0; }
gemma -bfile pruned_v2.1 -gk 1 -maf 0.05 -o kin
