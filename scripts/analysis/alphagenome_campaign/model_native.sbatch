#!/bin/bash
#SBATCH --job-name=masld_ag_native_full
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_native_full_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_native_full_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba
export XLA_PYTHON_CLIENT_PREALLOCATE=false HF_HUB_OFFLINE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python
OUT=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model/native
mkdir -p "$OUT/env"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv > "$OUT/env/gpu_${SLURM_JOB_ID}.csv"
"$PY" -m pip freeze > "$OUT/env/packages_${SLURM_JOB_ID}.txt"
if [[ ! -e "$OUT/currin_variants.tsv" ]]; then
  "$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_baseline.py" prepare --out "$OUT"
fi
"$PY" "$PROJ/scripts/analysis/alphagenome_program/i1_01_score_local.py" \
  --variants "$OUT/currin_variants.tsv" --length 1048576 \
  --out "$OUT/native_1048576.tsv" --resume --require-card l40s --flush-every 20
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_baseline.py" evaluate --out "$OUT"
