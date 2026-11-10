#!/bin/bash
# Process ONE accession end to end:
#   stream unplaced-unmapped reads from the public BAM -> remove organelle/PhiX
#   -> Kraken2 + Bracken -> per-accession QC line. Deletes FASTQ when done.
# Usage: scripts/01_process_accession.sh <ACCESSION> <THREADS> [--fetch-only|--classify-only]
set -euo pipefail
ACC=$1; T=${2:-8}; MODE=${3:-all}
source "$(dirname "$0")/../config.sh"
R="$WORK/results"; FQ="$WORK/fastq/$ACC"; mkdir -p "$FQ"
cd "$FQ"   # samtools caches the remote .bai in CWD; keep it in the per-accession dir (deleted at end)
[ -s "$R/qc/$ACC.qc.tsv" ] && { echo "$ACC already done"; exit 0; }
# Some deposited BAMs have no unplaced-unmapped block (unmapped pairs stripped): skip, don't retry.
NUNP=$(awk -F'\t' -v a="$ACC" '$1==a{print $4}' "$R/qc/total_reads.tsv" 2>/dev/null)
if [ "${NUNP:-1}" = "0" ]; then
  grep -qx "$ACC" "$PROJ/logs/no_unmapped.txt" 2>/dev/null || echo "$ACC" >> "$PROJ/logs/no_unmapped.txt"
  echo "$ACC has no unmapped block in its BAM; skipped"; rm -rf "$FQ"; exit 0
fi

# ---- 1. fetch: only the '*' (unplaced unmapped) block, via the BAM index ----
if [ "$MODE" != "--classify-only" ] && [ ! -e "$FQ/done.fetch" ]; then
  for try in 1 2 3; do
    if samtools view -u "$BAM_URL/$ACC.realigned.bam" '*' \
       | samtools fastq -n -1 "$FQ/r1.fq.gz" -2 "$FQ/r2.fq.gz" -0 /dev/null -s /dev/null - ; then
      touch "$FQ/done.fetch"; break; fi
    echo "fetch failed ($try), retrying"; sleep $((try*30))
  done
  [ -e "$FQ/done.fetch" ] || { echo "FETCH FAILED $ACC" >&2; exit 1; }
fi
[ "$MODE" = "--fetch-only" ] && exit 0

# ---- 2. decoy removal: rice chloroplast + mitochondrion + PhiX ----
# Keep reads hitting the mitochondrion as BAM (input for the CMS / cytoplasm follow-up).
minimap2 -t "$T" -ax sr "$WORK/refs/decoy.fa" "$FQ/r1.fq.gz" "$FQ/r2.fq.gz" 2>/dev/null \
  | samtools view -@2 -b -o "$FQ/decoy.bam" -
samtools sort -@2 -o "$FQ/decoy.sorted.bam" "$FQ/decoy.bam" && samtools index "$FQ/decoy.sorted.bam"
samtools view -b -F 0x904 "$FQ/decoy.sorted.bam" BA000029.3 > "$R/mito/$ACC.mito.bam"
read -r N_CP N_MT N_PHIX < <(samtools idxstats "$FQ/decoy.sorted.bam" \
  | awk '$1~/NC_001320/{c=$3} $1~/BA000029/{m=$3} $1~/NC_001422/{p=$3} END{print c+0,m+0,p+0}')
# pairs where neither mate hit a decoy (flag 12 = both unmapped)
samtools fastq -n -f 12 -F 0x900 -1 "$FQ/c1.fq" -2 "$FQ/c2.fq" -0 /dev/null -s /dev/null "$FQ/decoy.bam" 2>/dev/null

# ---- 3. Kraken2 + Bracken ----
kraken2 --db "$KDB" --threads "$T" --paired --confidence "$K2_CONF" \
        --report-minimizer-data --report "$R/kraken/$ACC.k2report" \
        "$FQ/c1.fq" "$FQ/c2.fq" > /dev/null 2> "$R/kraken/$ACC.k2log"
READLEN=$(seqkit head -n 20000 "$FQ/c1.fq" | seqkit fx2tab -n -l | awk '{s+=$NF} END{l=int(s/NR); print (l>=140?150:(l>=90?100:75))}')
for LV in S G; do
  bracken -d "$KDB" -i "$R/kraken/$ACC.k2report" -o "$R/bracken/$ACC.$LV.bracken" \
          -r "$READLEN" -l "$LV" -t 5 > /dev/null 2>&1 || echo "bracken $LV failed $ACC" >&2
done

# ---- 3b. independent check: direct alignment to pathogen genomes (breadth of coverage) ----
minimap2 -t "$T" -ax sr "$WORK/refs/pathogens.mmi" "$FQ/c1.fq" "$FQ/c2.fq" 2>/dev/null \
  | samtools sort -@2 -o "$FQ/patho.bam" - && samtools index "$FQ/patho.bam"
samtools coverage --min-MQ 20 "$FQ/patho.bam" | sed "s/^/$ACC\t/" | tail -n +2 > "$R/patho/$ACC.coverage.tsv"

# ---- 4. QC line ----
N_UNM=$(( $(zcat "$FQ/r1.fq.gz" | wc -l) / 4 * 2 ))
N_CLEAN=$(( $(wc -l < "$FQ/c1.fq") / 4 * 2 ))
N_CLASS=$(awk -F'\t' '$7==1{print $2}' "$R/kraken/$ACC.k2report")   # root clade (minimizer-data report: taxid = col 7)
GC=$(seqkit fx2tab -n -g "$FQ/c1.fq" | awk '{s+=$NF} END{printf "%.1f", s/NR}')
echo -e "$ACC\t$N_UNM\t$N_CP\t$N_MT\t$N_PHIX\t$N_CLEAN\t${N_CLASS:-0}\t$READLEN\t$GC" > "$R/qc/$ACC.qc.tsv"

rm -rf "$FQ"
echo "$ACC done"
