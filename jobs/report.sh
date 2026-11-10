#!/bin/bash
# Scheduler-independent job body; submitted by run_pipeline.sh (PBS, SLURM or local).
# Receives PROJ, WORK, LIST, NCPU (and TASK where noted) from the environment.
set -uo pipefail
cd "$PROJ" && source config.sh && activate_env || exit 1
mkdir -p "$PROJ/logs"
# Stage 7: figures, GWAS summary, locus robustness check, and a copy of the small results into the repo.
set -e
if ls "$WORK"/gwas/output/lmm_*.assoc.txt >/dev/null 2>&1; then
  python scripts/04_plots.py --gwas
  python scripts/05_locus_check.py present_Xanthomonas_oryzae || echo "locus check skipped"
else
  python scripts/04_plots.py
fi
bash scripts/06_collect_results.sh
