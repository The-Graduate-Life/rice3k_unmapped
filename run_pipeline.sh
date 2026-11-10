#!/bin/bash
# ---------------------------------------------------------------------------
# Run the whole analysis with one command. Every stage is submitted at once and chained with job
# dependencies, so you can log out after this returns.
#
#   bash run_pipeline.sh                       # full run, scheduler auto-detected (PBS / SLURM / local)
#   bash run_pipeline.sh --from gwas           # resume from a stage
#   bash run_pipeline.sh --status              # progress report
#   bash run_pipeline.sh --dry-run             # print the submissions, change nothing
#
# Stages:  setup -> process (N jobs) -> aggregate -> identity (N jobs) -> aggregate
#                -> kinship + gwas (N jobs) -> report
# All stages skip work that is already done, so re-running is always safe.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")"
source config.sh

usage() { sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; cat <<EOF

Options:
  --scheduler pbs|slurm|local   force a scheduler (default: auto-detect)
  --list FILE                   accession list (default: metadata/accessions_all.txt)
  --from STAGE                  start at setup|process|identity|gwas|report (default: setup)
  --no-gwas                     stop after the presence calls (no kinship/GWAS/locus check)
  --classify-only               reads were pre-fetched on a login node (no internet on compute nodes)
  --check-node                  only submit a short job that tests a compute node
  --status                      show progress and exit
  --dry-run                     print what would be submitted
EOF
}

SCHED=auto; LIST="$PROJ/metadata/accessions_all.txt"; FROM=setup; GWAS=1; DRY=0; CHECK=0; STATUS=0; MODE=all
while [ $# -gt 0 ]; do
  case $1 in
    --scheduler) SCHED=$2; shift ;;
    --list) LIST=$(readlink -f "$2"); shift ;;
    --from) FROM=$2; shift ;;
    --no-gwas) GWAS=0 ;;
    --classify-only) MODE=--classify-only ;;
    --check-node) CHECK=1 ;;
    --status) STATUS=1 ;;
    --dry-run) DRY=1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown option $1"; usage; exit 1 ;;
  esac; shift
done
export LIST
LOGS="$PROJ/logs"; mkdir -p "$LOGS"

# ---------------------------------------------------------------- status
if [ $STATUS = 1 ]; then
  N=$(grep -c . "$LIST")
  q() { ls $1 2>/dev/null | wc -l; }
  echo "accessions in list        : $N"
  echo "processed (qc files)      : $(q "$WORK/results/qc/*.qc.tsv")"
  echo "skipped, no unmapped block: $( [ -s "$LOGS/no_unmapped.txt" ] && wc -l < "$LOGS/no_unmapped.txt" || echo 0)"
  echo "failed (logs/failed.txt)  : $( [ -s "$LOGS/failed.txt" ] && sort -u "$LOGS/failed.txt" | wc -l || echo 0)"
  echo "identity checks done      : $(q "$WORK/results/identity/*.tsv")"
  echo "GWAS scans done           : $(q "$WORK/gwas/output/lmm_*.assoc.txt")"
  echo "figures                   : $(q "$WORK/results/figures/*.png")"
  command -v qstat >/dev/null 2>&1 && { echo; qstat -u "$USER" 2>/dev/null | tail -n +4 || true; }
  command -v squeue >/dev/null 2>&1 && { echo; squeue -u "$USER" 2>/dev/null || true; }
  exit 0
fi

# ---------------------------------------------------------------- scheduler
if [ "$SCHED" = auto ]; then
  if command -v qsub >/dev/null 2>&1 && qstat --version 2>&1 | grep -qi pbs; then SCHED=pbs
  elif command -v sbatch >/dev/null 2>&1; then SCHED=slurm
  elif command -v qsub >/dev/null 2>&1; then SCHED=pbs
  else SCHED=local; fi
fi
echo "scheduler: $SCHED | list: $LIST ($(grep -c . "$LIST") accessions) | WORK: $WORK | Kraken2 DB: $KDB_FLAVOR"

# submit NAME RESOURCES "DEP_IDS" SCRIPT [VAR=value ...]  -> prints the job id
submit() {
  local name=$1 res=$2 deps=$3 script=$4; shift 4
  local q cpu mem wall; IFS='|' read -r q cpu mem wall <<< "$res"
  local vars="PROJ=$PROJ,WORK=$WORK,LIST=$LIST,NCPU=$cpu"
  local kv; for kv in "$@"; do vars+=",$kv"; done
  local cmd
  case $SCHED in
    pbs)   cmd=(qsub -N "$name" ${q:+-q "$q"} -l "nodes=1:ppn=$cpu" -l "mem=${mem}gb" -l "walltime=$wall"
                -j oe -o "$LOGS/$name.log" -v "$vars" ${deps:+-W "depend=afterok:$deps"} "$script") ;;
    slurm) cmd=(sbatch --parsable -J "$name" ${q:+-p "$q"} -c "$cpu" --mem="${mem}G" -t "$wall"
                -o "$LOGS/$name.log" --export="ALL,$vars" ${deps:+--dependency="afterok:$deps"} "$script") ;;
    local) cmd=(env NCPU="$(nproc)" "$@" bash "$script") ;;
  esac
  if [ $DRY = 1 ]; then echo "  [dry-run] ${cmd[*]}" >&2; echo "JOB_$name"; return; fi
  if [ $SCHED = local ]; then
    echo "  running $name ..." >&2
    "${cmd[@]}" > "$LOGS/$name.log" 2>&1 || { echo "  $name FAILED: see $LOGS/$name.log" >&2; exit 1; }
    echo "local"
  else
    "${cmd[@]}"
  fi
}
join() { local IFS=:; echo "$*" | sed 's/local//g; s/::*/:/g; s/^://; s/:$//'; }

if [ $CHECK = 1 ]; then
  id=$(submit rice3k_check "$RES_CHECK" "" jobs/check_node.sh)
  echo "check job $id -> $LOGS/rice3k_check.log"; exit 0
fi

stage_n() { case $1 in setup) echo 1;; process) echo 2;; identity) echo 3;; gwas) echo 4;; report) echo 5;; *) echo "bad --from $1" >&2; exit 1;; esac; }
S=$(stage_n "$FROM")
DEP=""

# 1. setup ---------------------------------------------------------------------
if [ "$S" -le 1 ]; then
  DEP=$(submit rice3k_setup "$RES_SETUP" "" jobs/setup.sh); echo "setup     : $DEP"
fi
SETUP_DEP=$DEP

# 2. process + aggregate ----------------------------------------------------------
if [ "$S" -le 2 ]; then
  NT=$(( ( $(grep -c . "$LIST") + BATCH - 1 ) / BATCH )); IDS=()
  for i in $(seq 1 "$NT"); do IDS+=("$(submit "rice3k_proc_$i" "$RES_PROCESS" "$DEP" jobs/process.sh TASK="$i" MODE="$MODE")"); done
  echo "process   : $NT jobs"
  DEP=$(submit rice3k_agg1 "$RES_AGGREGATE" "$(join "${IDS[@]}")" jobs/aggregate.sh); echo "aggregate : $DEP"
fi

# 3. identity + aggregate ----------------------------------------------------------
if [ "$S" -le 3 ]; then
  IDS=()
  for i in $(seq 1 "$N_IDENTITY_JOBS"); do IDS+=("$(submit "rice3k_id_$i" "$RES_IDENTITY" "$DEP" jobs/identity.sh TASK="$i")"); done
  echo "identity  : $N_IDENTITY_JOBS jobs"
  DEP=$(submit rice3k_agg2 "$RES_AGGREGATE" "$(join "${IDS[@]}")" jobs/aggregate.sh); echo "aggregate : $DEP"
fi

# 4. kinship + GWAS ------------------------------------------------------------------
if [ $GWAS = 1 ] && [ "$S" -le 4 ]; then
  KJ=$(submit rice3k_kin "$RES_KINSHIP" "$SETUP_DEP" jobs/kinship.sh); echo "kinship   : $KJ"
  IDS=()
  for i in $(seq 1 "$N_GWAS_JOBS"); do IDS+=("$(submit "rice3k_gwas_$i" "$RES_GWAS" "$(join "$DEP" "$KJ")" jobs/gwas.sh TASK="$i")"); done
  echo "gwas      : $N_GWAS_JOBS jobs"
  DEP=$(join "${IDS[@]}")
fi

# 5. report ------------------------------------------------------------------------
REP=$(submit rice3k_report "$RES_REPORT" "$DEP" jobs/report.sh); echo "report    : $REP"

echo
if [ $SCHED = local ]; then
  echo "Finished. Results: $WORK/results and ${RESULTS_DIR:-$PROJ/results}"
else
  echo "All jobs submitted. Check progress any time with:  bash run_pipeline.sh --status"
  echo "Logs: $LOGS/<job name>.log (PBS writes them when each job ends). Final results land in $PROJ/results/."
fi
