#!/bin/bash
# Scheduler-independent job body; submitted by run_pipeline.sh (PBS, SLURM or local).
# Receives PROJ, WORK, LIST, NCPU (and TASK where noted) from the environment.
set -uo pipefail
cd "$PROJ" && source config.sh && activate_env || exit 1
mkdir -p "$PROJ/logs"
# Probe a compute node: RAM, scratch, and whether it can reach the public S3 bucket.
echo "node: $(hostname)"; free -g | head -2; mkdir -p "$WORK"; df -h "$WORK"
if samtools view "$BAM_URL/B001.realigned.bam" '*' 2>/dev/null | head -1 | grep -q .; then
  echo "INTERNET_FROM_COMPUTE=yes"
else
  echo "INTERNET_FROM_COMPUTE=no  -> see docs/HPC_GUIDE.md (fetch on a login node first)"
fi
rm -f B001.realigned.bam.bai
