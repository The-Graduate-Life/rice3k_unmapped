#!/bin/bash
# One-command software setup: installs Miniforge (if no conda/mamba is found) and creates the
# "reuse3k" conda environment from env/reuse3k.yml. Safe to re-run.
#   bash install.sh
set -euo pipefail
cd "$(dirname "$0")"
ENV_NAME=reuse3k

for c in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3"; do
  [ -s "$c/etc/profile.d/conda.sh" ] && { source "$c/etc/profile.d/conda.sh"; break; }
done

if ! command -v conda >/dev/null 2>&1; then
  echo ">> no conda found: installing Miniforge into $HOME/miniforge3"
  OS=$(uname -s); ARCH=$(uname -m)
  curl -fsSL -o /tmp/miniforge.sh \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-${OS}-${ARCH}.sh"
  bash /tmp/miniforge.sh -b -p "$HOME/miniforge3" && rm /tmp/miniforge.sh
  source "$HOME/miniforge3/etc/profile.d/conda.sh"
fi

SOLVER=conda; command -v mamba >/dev/null 2>&1 && SOLVER=mamba
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo ">> updating existing environment $ENV_NAME"
  $SOLVER env update -n "$ENV_NAME" -f env/reuse3k.yml --prune
else
  echo ">> creating environment $ENV_NAME"
  $SOLVER env create -n "$ENV_NAME" -f env/reuse3k.yml
fi

conda activate "$ENV_NAME"
echo ">> tool check"
for t in samtools minimap2 seqkit kraken2 bracken plink gemma python curl; do
  printf '   %-9s %s\n' "$t" "$(command -v $t >/dev/null && echo ok || echo MISSING)"
done
echo
echo "Done. Next: bash run_pipeline.sh --help   (or follow FOLLOWME.md)"
