#!/bin/bash
# Copy the small, shareable results from $WORK into the repository (results/), plus a record of the
# exact software versions. Large intermediates (FASTQ, Kraken reports, GEMMA association files)
# stay in $WORK.
#   bash scripts/06_collect_results.sh
set -euo pipefail
source "$(dirname "$0")/../config.sh"
R="$WORK/results"; D="${RESULTS_DIR:-$PROJ/results}"
mkdir -p "$D/tables" "$D/figures" "$D/gwas"

for f in qc_summary.tsv targets.tsv target_hits.tsv prevalence_by_subpop.tsv prevalence_by_country.tsv \
         prevalence_by_batch.tsv gwas_summary.tsv gwas_pheno_names.txt; do
  [ -s "$R/tables/$f" ] && cp "$R/tables/$f" "$D/tables/"
done
cp "$R"/tables/gwas_hits_*.tsv "$R"/tables/locus_check_*.tsv "$D/gwas/" 2>/dev/null || true
cp "$R"/figures/*.png "$D/figures/" 2>/dev/null || true
for f in no_unmapped.txt failed.txt identity_failed.txt; do
  [ -s "$PROJ/logs/$f" ] && cp "$PROJ/logs/$f" "$D/tables/run_$f"
done
# species x accession matrix is ~10-30 MB: compress it
[ -s "$R/tables/species_rpm.tsv" ] && gzip -c "$R/tables/species_rpm.tsv" > "$D/tables/species_rpm.tsv.gz"

# provenance
if command -v conda >/dev/null 2>&1; then conda env export --no-builds 2>/dev/null | grep -v '^prefix:' > "$D/reuse3k.lock.yml" || true; fi
{
  echo "collected: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host: $(hostname)"
  echo "kraken2_db: $KDB_FLAVOR"
  echo "accessions_with_results: $(ls "$R"/qc/*.qc.tsv 2>/dev/null | wc -l)"
  v() { printf '%s: %s\n' "$1" "$(eval "$2" 2>&1 | grep -m1 -iE "$3" | tr -d '\r' || true)"; }
  v samtools "samtools --version" samtools
  v minimap2 "minimap2 --version" .
  v kraken2  "kraken2 --version" version
  v bracken  "bracken -v" bracken
  v seqkit   "seqkit version" seqkit
  v gemma    "gemma -h" version
  v plink    "plink --version" plink
  v python   "python --version" python
} > "$D/RUN_INFO.txt"
echo "results copied to $D (see $D/RUN_INFO.txt)"
