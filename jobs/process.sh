#!/bin/bash
# Scheduler-independent job body; submitted by run_pipeline.sh (PBS, SLURM or local).
# Receives PROJ, WORK, LIST, NCPU (and TASK where noted) from the environment.
set -uo pipefail
cd "$PROJ" && source config.sh && activate_env || exit 1
mkdir -p "$PROJ/logs"
# Stage 2: TASK = batch number; processes lines (TASK-1)*BATCH+1 .. TASK*BATCH of LIST.
# A failed accession is logged to logs/failed.txt and does not stop the batch.
START=$(( (TASK - 1) * BATCH + 1 )); END=$(( START + BATCH - 1 ))
echo "process task $TASK: lines $START-$END of $LIST on $(hostname), $NCPU cpus"
sed -n "${START},${END}p" "$LIST" | while read -r ACC; do
  [ -z "$ACC" ] && continue
  bash scripts/01_process_accession.sh "$ACC" "$NCPU" "${MODE:-all}" \
    || echo -e "$ACC\tFAILED" >> "$PROJ/logs/failed.txt"
done
echo "process task $TASK finished $(date)"
