#!/bin/bash
#SBATCH --job-name=masld_ag_scalar_pilot
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_scalar_pilot_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_scalar_pilot_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba
export XLA_PYTHON_CLIENT_PREALLOCATE=false HF_HUB_OFFLINE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model
CONFIG=${MODEL_CONFIG:-$BASE/configs/adapter_r4_last3_variant.json}
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_scalar_tests.py" --out "$BASE/inspection/scalar_invariants_gpu_${SLURM_JOB_ID}.json"
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_train.py" --config "$CONFIG" \
  --out "$BASE/pilot_${SLURM_JOB_ID}" --steps 100 --microbatch 1 --validation-limit 32 --max-hours 3.8
