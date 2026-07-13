#!/bin/bash
#SBATCH --job-name=scvijack
#SBATCH --partition=gpu
#SBATCH --qos=interactive
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-48
#SBATCH --output=Analysis/SingleCell/scripts/logs_scvi_loocv/343v_arr_%A_%a.log
#SBATCH --error=Analysis/SingleCell/scripts/logs_scvi_loocv/343v_arr_%A_%a.err
#
# Leakage-free per-donor scVI F-stage jackknife (audit P0#1, 2026-06-29).
# Retrains scVI FROM SCRATCH per held-out donor (latent-only projection of the
# held-out cells) -> kNN-5 on documented Andrews F-stage. NOT the fixed-latent
# head (that leak gives the inflated QWK 0.74-0.76). Aggregate with 343x.
#
# Array maps SLURM_ARRAY_TASK_ID -> line N of 343v_remaining_donors.txt (the 48
# documented-F-stage donors not already on disk; 10 of the 58 were done prior).

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/logs_scvi_loocv

DONOR_LIST=Analysis/SingleCell/scripts/343v_remaining_donors.txt
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
