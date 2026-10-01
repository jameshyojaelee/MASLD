#!/bin/bash
#SBATCH --job-name=masld_ag_frozen_heads
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --dependency=afterok:21772543
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_frozen_heads_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_frozen_heads_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_frozen_fit.py" --out "$BASE/frozen" --configs "$BASE/configs"
