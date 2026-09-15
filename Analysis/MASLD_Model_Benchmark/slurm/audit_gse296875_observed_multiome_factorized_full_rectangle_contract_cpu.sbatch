#!/bin/bash
#SBATCH --job-name=model-check-321
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON="$ROOT/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
cd "$ROOT"
mkdir -p executions

PYTHONPATH="src:." "$PYTHON" -m unittest -v \
  tests.unit.test_observed_multiome_factorized_full_rectangle_contract
PYTHONPATH="src:." "$PYTHON" scripts/audit_gse296875_observed_multiome_factorized_full_rectangle_contract.py \
  --root "$ROOT" \
  --config "$ROOT/config/campaigns/gse296875_observed_multiome_factorized_full_rectangle_contract_20260825.json" \
  --output "$ROOT/executions/model-check-321-${SLURM_JOB_ID}"
