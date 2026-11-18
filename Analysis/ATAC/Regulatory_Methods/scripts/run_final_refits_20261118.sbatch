#!/usr/bin/env bash
#SBATCH --job-name=ws4a_refit_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --array=0-4
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=Analysis/ATAC/Regulatory_Methods/logs/ws4a_refit_20261118_%A_%a.out
#SBATCH --error=Analysis/ATAC/Regulatory_Methods/logs/ws4a_refit_20261118_%A_%a.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
PY_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python"

"$PY_BIN" Analysis/ATAC/Regulatory_Methods/scripts/final_prediction_refits_20261118.py --workers 8
