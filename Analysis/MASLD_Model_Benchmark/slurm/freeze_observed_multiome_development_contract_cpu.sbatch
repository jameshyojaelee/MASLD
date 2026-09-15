#!/bin/bash
#SBATCH --job-name=model-check-255
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
cd "$ROOT"
mkdir -p logs executions

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_observed_multiome_development_contract
PYTHONPATH="src:." python scripts/freeze_observed_multiome_development_contract.py \
  --root "$ROOT" \
  --output "$ROOT/executions/model-check-255-${SLURM_JOB_ID}"
