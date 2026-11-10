# shellcheck shell=bash
# ---------------------------------------------------------------------------
# rice3k_unmapped configuration. Sourced by every script and job.
# Every value can be overridden from the environment, e.g.
#     WORK=/big/disk/rice3k_work KDB_FLAVOR=viral bash run_pipeline.sh ...
# ---------------------------------------------------------------------------

# ---- paths ----------------------------------------------------------------
# PROJ = this repository (code + metadata). WORK = large scratch area for downloads and results
# (~45 GB peak with the default Kraken2 database).
export PROJ="${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
export WORK="${WORK:-${SCRATCH:-$HOME}/rice3k_work}"

# ---- software environment -------------------------------------------------
# install.sh creates a conda environment called "reuse3k". If yours lives elsewhere, set CONDA_ENV
# to its full path (e.g. /home/me/miniforge3/envs/reuse3k).
export CONDA_ENV="${CONDA_ENV:-reuse3k}"

# ---- analysis parameters --------------------------------------------------
# Kraken2 database: pluspf16 (12 GB, ~20 GB RAM; used for the published run) | pluspf (83 GB, ~120 GB RAM)
#                   | viral (0.5 GB; smoke test only)
export KDB_FLAVOR="${KDB_FLAVOR:-pluspf16}"
export KDB="$WORK/db/kraken_${KDB_FLAVOR}"
export K2_CONF="${K2_CONF:-0.1}"           # Kraken2 --confidence
export BATCH="${BATCH:-25}"               # accessions per processing job
export N_IDENTITY_JOBS="${N_IDENTITY_JOBS:-6}"
export N_GWAS_JOBS="${N_GWAS_JOBS:-12}"   # traits are spread over this many jobs

# ---- data sources (public, no account needed) -----------------------------
export BAM_URL="https://3kricegenome.s3.us-east-1.amazonaws.com/Nipponbare"
export S3_BASE="https://3kricegenome.s3.amazonaws.com"

# ---- cluster resources: "queue|cpus|mem_GB|walltime" ------------------------
# Defaults are for the Alabama Supercomputer (PBS Pro: small <=4 GB/8 cpu, medium <=16 GB/16 cpu,
# large <=120 GB/128 cpu, every queue >= 1 h). On SLURM, "queue" is the partition; leave it empty
# to use the cluster default.
export RES_CHECK="${RES_CHECK:-small|2|4|01:00:00}"
export RES_SETUP="${RES_SETUP:-medium|2|8|06:00:00}"
export RES_PROCESS="${RES_PROCESS:-large|8|24|08:00:00}"
export RES_AGGREGATE="${RES_AGGREGATE:-medium|2|16|02:00:00}"
export RES_IDENTITY="${RES_IDENTITY:-large|8|16|08:00:00}"
export RES_KINSHIP="${RES_KINSHIP:-medium|4|16|06:00:00}"
export RES_GWAS="${RES_GWAS:-medium|4|16|24:00:00}"
export RES_REPORT="${RES_REPORT:-medium|2|16|02:00:00}"

# ---- helpers ----------------------------------------------------------------
activate_env() {
  # already inside the environment?
  command -v kraken2 >/dev/null 2>&1 && command -v gemma >/dev/null 2>&1 && return 0
  if ! command -v conda >/dev/null 2>&1; then
    for c in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3"; do
      [ -s "$c/etc/profile.d/conda.sh" ] && { source "$c/etc/profile.d/conda.sh"; break; }
    done
  fi
  command -v conda >/dev/null 2>&1 || { echo "conda not found: run 'bash install.sh' first" >&2; return 1; }
  [ -n "$(type -t conda)" ] && [ "$(type -t conda)" != function ] && eval "$(conda shell.bash hook)"
  conda activate "$CONDA_ENV"
}
