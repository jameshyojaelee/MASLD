#!/bin/bash
#SBATCH --job-name=masld_inference_cost
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_inference_cost_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_inference_cost_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false HF_HUB_OFFLINE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=$PROJ/GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/inference_cost_${SLURM_JOB_ID}
/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python \
  "$PROJ/scripts/analysis/masld_two_models/inference_cost.py" --out "$OUT"
