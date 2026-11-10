#!/bin/bash
# Scheduler-independent job body; submitted by run_pipeline.sh (PBS, SLURM or local).
# Receives PROJ, WORK, LIST, NCPU (and TASK where noted) from the environment.
set -uo pipefail
cd "$PROJ" && source config.sh && activate_env || exit 1
mkdir -p "$PROJ/logs"
# Stage 6: GEMMA linear mixed model. TASK = 1..N_GWAS_JOBS; each job runs every N-th trait,
# once with covariates PC1-5 + depth ("base") and once adding sample-batch dummies ("batch").
# -lmm 1 = Wald test (p_wald), the statistic used downstream.
T="$WORK/results/tables"; cd "$WORK/gwas" || exit 1
while read -r IDX NAME; do
  [ $(( (IDX - 1) % N_GWAS_JOBS )) -eq $(( TASK - 1 )) ] || continue
  NN=$(cut -f"$IDX" "$T/gwas_pheno.txt" | grep -vc NA || true)
  [ "$NN" -lt 200 ] && { echo "$NAME: only $NN non-missing values, skipped"; continue; }
  for COVSET in base batch; do
    if [ $COVSET = batch ]; then COV="$T/gwas_covar_batch.txt"; SUF=_batch; else COV="$T/gwas_covar.txt"; SUF=; fi
    OUT="lmm_${NAME}${SUF}"
    [ -s "output/$OUT.assoc.txt" ] && { echo "$OUT exists"; continue; }
    echo "== $OUT"
    gemma -bfile pruned_v2.1 -k output/kin.cXX.txt -p "$T/gwas_pheno.txt" -n "$IDX" \
          -c "$COV" -lmm 1 -maf 0.05 -miss 0.1 -o "$OUT" > "output/$OUT.gemma.log" 2>&1 \
      || echo "GEMMA failed for $OUT (see $WORK/gwas/output/$OUT.gemma.log)"
  done
done < "$T/gwas_pheno_names.txt"
echo "gwas task $TASK finished $(date)"
