#!/bin/bash
#SBATCH --job-name=343w_coh
#SBATCH --partition=gpu
#SBATCH --qos=interactive
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=Analysis/SingleCell/scripts/logs_scvi_loocv/343w_%x_%j.log
#SBATCH --error=Analysis/SingleCell/scripts/logs_scvi_loocv/343w_%x_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation

eval "$(micromamba shell hook --shell bash)"
micromamba activate rapids_singlecell

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "[run] cohort=${COHORT:?must export COHORT=<dataset_name>}"
nvidia-smi || true

python Analysis/SingleCell/scripts/343w_scvi_cohort_out.py --hold-out-cohort "${COHORT}"
