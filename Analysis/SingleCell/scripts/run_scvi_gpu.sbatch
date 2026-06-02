#!/usr/bin/env bash
#SBATCH --job-name=scvi_int_gpu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/logs/scvi_gpu_%j.out
#SBATCH --error=Analysis/SingleCell/logs/scvi_gpu_%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CONDA_ENV="${PROJECT_ROOT}/.agent/scanpy_env"
SCRIPT="${PROJECT_ROOT}/Analysis/SingleCell/scripts/scvi_integration.py"

mkdir -p "${PROJECT_ROOT}/Analysis/SingleCell/logs"

eval "$(micromamba shell hook -s bash)"
micromamba activate "${CONDA_ENV}"

# Fix: conda env's libstdc++ must take precedence over system/module versions
export LD_LIBRARY_PATH="${CONDA_ENV}/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1

# Pre-download CellTypist models if not cached
python3 -c "import celltypist; celltypist.models.download_models(force_update=False)" 2>/dev/null || true

# Unbuffered output so SLURM logs update in real-time
export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1

# Run with isolated output directory
python3 -X faulthandler "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --output-dir "Analysis/SingleCell/results_gpu" \
    --export-data "${PROJECT_ROOT}/Analysis/SingleCell/results_gpu/pre_scvi_checkpoint.h5ad"
