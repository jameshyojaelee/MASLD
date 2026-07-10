#!/bin/bash
#SBATCH --job-name=scvijack
#SBATCH --partition=gpu
#SBATCH --qos=interactive
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-10
#SBATCH --output=Analysis/SingleCell/scripts/logs_scvi_loocv/343v_stale10_%A_%a.log
#SBATCH --error=Analysis/SingleCell/scripts/logs_scvi_loocv/343v_stale10_%A_%a.err
#
# Regenerate the 10 ORIGINAL (May-20, pre-current-script) donor-jackknife folds
# with the current leakage-free 343v so the full 58-donor QWK is single-method.
# (Audit verification 2026-06-30: the 48 "remaining" folds were refreshed 06-29;
#  these 10 were left over from the original run. Refresh for a clean aggregate.)

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/logs_scvi_loocv

DONOR_LIST=Analysis/SingleCell/scripts/343v_stale10_donors.txt
DONOR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$DONOR_LIST")
if [ -z "${DONOR:-}" ]; then
  echo "[fatal] no donor at line ${SLURM_ARRAY_TASK_ID} of $DONOR_LIST" >&2
  exit 1
fi

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rapids_singlecell

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "[run] task=${SLURM_ARRAY_TASK_ID} donor=${DONOR}"
nvidia-smi || true

python Analysis/SingleCell/scripts/343v_scvi_donor_jackknife.py --hold-out-donor "${DONOR}"
