#!/bin/bash
# Scheduler-independent job body; submitted by run_pipeline.sh (PBS, SLURM or local).
# Receives PROJ, WORK, LIST, NCPU (and TASK where noted) from the environment.
set -uo pipefail
cd "$PROJ" && source config.sh && activate_env || exit 1
mkdir -p "$PROJ/logs"
# Stage 4: per-read identity for every accession with a "confirmed" or "congener" call.
# TASK = 1..N_IDENTITY_JOBS; each job takes every N-th accession of the list.
set -e
HITS="$WORK/results/tables/target_hits.tsv"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="call") k=i; next} $k=="confirmed" || $k=="congener" {print $1}' "$HITS" \
  | sort -u | awk -v n="$N_IDENTITY_JOBS" -v t="$TASK" '(NR-1) % n == t-1' > "$WORK/identity_list_$TASK.txt"
echo "identity task $TASK: $(wc -l < "$WORK/identity_list_$TASK.txt") accessions"
set +e
while read -r ACC; do
  bash scripts/02_identity_check.sh "$ACC" "$NCPU" || echo "$ACC" >> "$PROJ/logs/identity_failed.txt"
done < "$WORK/identity_list_$TASK.txt"
echo "identity task $TASK finished $(date)"
