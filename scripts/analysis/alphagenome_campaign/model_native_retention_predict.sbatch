#!/bin/bash
#SBATCH --job-name=masld_native_ATAC_retention
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --dependency=afterok:21773996
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_native_ATAC_retention_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_native_ATAC_retention_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export XLA_PYTHON_CLIENT_PREALLOCATE=false HF_HUB_OFFLINE=1
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
NATIVE_RETENTION_PANEL=$BASE/model/native_retention_panel_21773996
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_native_retention_predict.py" \
 --fixture "$NATIVE_RETENTION_PANEL" --sweep "$BASE/model/sweep_21772726" \
 --max-hours 1.95 --bootstrap 1000 --out "$BASE/model/native_retention_${SLURM_JOB_ID}"
