#!/bin/bash
# Per-read alignment identity of unmapped reads vs the pathogen panel.
# A true pathogen gives reads at >=97-99% identity; a related species (congener, endogenous
# viral sequence) aligns genome-wide but at ~80-92%. Breadth alone cannot tell these apart.
# Usage: scripts/02_identity_check.sh <ACCESSION> <THREADS>   -> results/identity/<ACC>.tsv
# columns: accession genome reads mean_identity frac_ge97 frac_ge99
set -euo pipefail
ACC=$1; T=${2:-8}
source "$(dirname "$0")/../config.sh"
OUT="$WORK/results/identity"; mkdir -p "$OUT"
[ -s "$OUT/$ACC.tsv" ] && { echo "$ACC already done"; exit 0; }
D="$WORK/idtmp/$ACC"; mkdir -p "$D"; cd "$D"
for try in 1 2 3; do
  if samtools view -u "$BAM_URL/$ACC.realigned.bam" '*' | samtools fastq -n - 2>/dev/null \
     | minimap2 -t "$T" -ax sr "$WORK/refs/pathogens.mmi" - 2>/dev/null \
     | samtools view -q 20 -F 0x904 - \
     | awk -F'\t' -v OFS='\t' -v acc="$ACC" '
         NR==FNR { g[$1]=$2; next }
         { c=$6; m=0
           while (match(c, /[0-9]+[MID]/)) { l=substr(c,RSTART,RLENGTH); m+=substr(l,1,length(l)-1); c=substr(c,RSTART+RLENGTH) }
           if (m < 50) next
           nm=0; for (i=12;i<=NF;i++) if ($i ~ /^NM:i:/) nm=substr($i,6)
           id=1-nm/m; G=g[$3]; n[G]++; s[G]+=id; if (id>=0.97) a[G]++; if (id>=0.99) b[G]++ }
         END { for (G in n) print acc, G, n[G], s[G]/n[G], a[G]/n[G], b[G]/n[G] }' \
       "$WORK/refs/pathogen_contigs.tsv" - > "$D/id.tsv"; then
    mv "$D/id.tsv" "$OUT/$ACC.tsv"; [ -s "$OUT/$ACC.tsv" ] || echo -e "$ACC\tnone\t0\t0\t0\t0" > "$OUT/$ACC.tsv"
    break
  fi
  echo "attempt $try failed"; sleep $((try*30))
done
cd "$WORK"; rm -rf "$D"
[ -s "$OUT/$ACC.tsv" ] || { echo "IDENTITY FAILED $ACC" >&2; exit 1; }
echo "$ACC done"
