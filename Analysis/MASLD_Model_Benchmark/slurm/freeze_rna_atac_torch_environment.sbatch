#!/usr/bin/env bash
#SBATCH --job-name=masld-rna-atac-env-freeze
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
PREFIX="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env"
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"

module load python/3.11.5-GCCcore-13.2.0
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "${ROOT}"

PYTHONPATH=src python scripts/freeze_python_environment.py \
  --micromamba "${MICROMAMBA}" \
  --prefix "${PREFIX}" \
  --runtime-id gpu_rna_atac_torch_smoke \
  --output-root executions/environment-locks
