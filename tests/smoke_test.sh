#!/bin/bash
# End-to-end plumbing test on 3 accessions, on any Linux/macOS machine (no cluster needed).
# Uses the small viral Kraken2 database (0.5 GB) and skips the SNP download and GWAS.
# Takes ~5-10 min and ~2 GB of disk; touches nothing in results/.
#   bash tests/smoke_test.sh            # WORK defaults to ./smoke_work (git-ignored)
set -euo pipefail
cd "$(dirname "$0")/.."
export WORK="${WORK:-$PWD/smoke_work}" KDB_FLAVOR=viral SKIP_GWAS_DATA=1 RESULTS_DIR="${WORK:-$PWD/smoke_work}/collected"
bash run_pipeline.sh --scheduler local --list tests/accessions_smoke.txt --no-gwas

echo; echo "== checks"
fail=0
check() { if eval "$2"; then echo "  ok   $1"; else echo "  FAIL $1"; fail=1; fi; }
T="$WORK/results/tables"
check "3 accessions processed"        '[ "$(ls "$WORK"/results/qc/*.qc.tsv | wc -l)" -eq 3 ]'
check "targets.tsv has 3 x 7 rows"    '[ "$(tail -n +2 "$T/targets.tsv" | wc -l)" -eq 21 ]'
check "calls are valid labels"        '! tail -n +2 "$T/targets.tsv" | awk -F"\t" "{print \$0}" | grep -qv -E "confirmed|congener|unvalidated|absent"'
check "5 QC/prevalence figures"       '[ "$(ls "$WORK"/results/figures/0*.png | wc -l)" -eq 5 ]'
check "mito BAMs written"             '[ "$(ls "$WORK"/results/mito/*.bam | wc -l)" -eq 3 ]'
check "results collected"             '[ -s "$RESULTS_DIR/RUN_INFO.txt" ]'
[ $fail = 0 ] && echo "SMOKE TEST PASSED" || { echo "SMOKE TEST FAILED (logs in logs/)"; exit 1; }
