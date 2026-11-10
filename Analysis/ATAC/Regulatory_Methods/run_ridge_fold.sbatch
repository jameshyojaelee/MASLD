#!/usr/bin/env bash
# RNA readout only; original outer/inner donor folds and feature selection.
#SBATCH --job-name=atac_ridge_grid
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods/logs/ridge_grid_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods/logs/ridge_grid_%j.log
set -euo pipefail
SOURCE=${1:?Frozen source directory required}
PREPARED=${2:?Verified prepared population required}
FOLD=${3:?Original outer fold required}
OUTPUT=${4:?New output required}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
MODEL_PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
/usr/bin/time -v "$MODEL_PY" "$SOURCE/ridge_grid_sensitivity.py" \
  --prepared "$PREPARED" --fold "$FOLD" --output "$OUTPUT" --workers "${SLURM_CPUS_PER_TASK:-8}"
