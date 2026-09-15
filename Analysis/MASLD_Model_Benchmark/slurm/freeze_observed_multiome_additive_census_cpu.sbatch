#!/bin/bash
#SBATCH --job-name=model-check-288
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=executions/%x-%j.out
#SBATCH --error=executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
cd "$ROOT"
mkdir -p executions

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_observed_multiome_additive_census
PYTHONPATH="src:." python scripts/freeze_observed_multiome_additive_census.py \
  --root "$ROOT" \
  --config "$ROOT/config/campaigns/observed_multiome_additive_census_revision_20260825.json" \
  --output "$ROOT/executions/model-check-288-${SLURM_JOB_ID}"
