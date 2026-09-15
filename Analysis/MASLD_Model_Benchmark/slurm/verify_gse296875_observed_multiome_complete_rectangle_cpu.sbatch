#!/usr/bin/env bash
#SBATCH --job-name=model-check-329
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON="$ROOT/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
cd "$ROOT"

PYTHONPATH="src:." "$PYTHON" -m unittest -v tests.unit.test_verify_observed_multiome_complete_rectangle tests.unit.test_stats
PYTHONPATH="src:." "$PYTHON" scripts/verify_gse296875_observed_multiome_complete_rectangle.py \
  --root "$ROOT" \
  --config "$ROOT/config/campaigns/gse296875_observed_multiome_complete_rectangle_verification_20260825.json" \
  --output "$ROOT/executions/model-check-329-${SLURM_JOB_ID}"
