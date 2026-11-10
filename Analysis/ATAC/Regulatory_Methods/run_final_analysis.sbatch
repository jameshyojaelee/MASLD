#!/usr/bin/env bash
#SBATCH --job-name=test
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods/logs/final_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods/logs/final_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
ANALYSIS_PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
if [[ "${1:-}" == "--command" ]]; then
  shift
  exec /usr/bin/time -v "$@"
fi
exec /usr/bin/time -v "$ANALYSIS_PY" "$@"
