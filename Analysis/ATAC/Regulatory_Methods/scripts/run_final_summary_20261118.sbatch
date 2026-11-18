#!/usr/bin/env bash
#SBATCH --job-name=ws4a_sum_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=Analysis/ATAC/Regulatory_Methods/logs/ws4a_sum_20261118_%j.out
#SBATCH --error=Analysis/ATAC/Regulatory_Methods/logs/ws4a_sum_20261118_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
PY_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python"

"$PY_BIN" Analysis/ATAC/Regulatory_Methods/scripts/final_prediction_summary_20261118.py
