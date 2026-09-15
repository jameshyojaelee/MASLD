#!/bin/bash
#SBATCH --job-name=model-cpu-train
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON="$ROOT/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
: "${SEED:?SEED must be exported by the frozen five-seed submission command}"
case "$SEED" in
  20260824|20260825|20260826|20260827|20260828) ;;
  *) printf 'Seed is outside frozen roster: %s\n' "$SEED" >&2; exit 2 ;;
esac
cd "$ROOT"
mkdir -p executions

PYTHONPATH="src:." "$PYTHON" -m unittest -v \
  tests.unit.test_observed_multiome_surface \
  tests.unit.test_observed_multiome_factorized \
  tests.unit.test_observed_multiome_factorized_seed_bundle
PYTHONPATH="src:." "$PYTHON" scripts/run_gse296875_observed_multiome_factorized_seed_bundle.py \
  --root "$ROOT" \
  --config "$ROOT/config/campaigns/gse296875_observed_multiome_factorized_full_rectangle_contract_20260825.json" \
  --contract-artifact "$ROOT/executions/model-check-321-21099900" \
  --seed "$SEED" \
  --output "$ROOT/executions/model-cpu-train-${SEED}-${SLURM_JOB_ID}"
