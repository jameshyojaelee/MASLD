#!/usr/bin/env bash
#SBATCH --job-name=model-training-069
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-training-068-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-training-068-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=3
export MKL_NUM_THREADS=3
export OPENBLAS_NUM_THREADS=3

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
CAMPAIGN=${ROOT}/executions/model-training-069-production-v2
MOLECULAR=${ROOT}/executions/model-data-064-21079902/fixture/molecular
FOLDS=${ROOT}/executions/model-data-064-21079902/fixture/folds
PREFLIGHT=${ROOT}/executions/model-training-067-21080776
DIAGNOSTIC=${ROOT}/executions/model-check-070-21082520
FIT_VIEWS=${PREFLIGHT}/fit_views
SURFACE=${ROOT}/config/evaluation/gse267145_histology_benchmark_surface.json
FITTER=${ROOT}/scripts/fit_gse267145_histology_baselines.py
AGGREGATOR=${ROOT}/scripts/aggregate_gse267145_histology_production_predictions.py
RUNNER=${ROOT}/scripts/run_gse267145_histology_production_bundle.py

cd "${ROOT}"
PYTHONPATH=src "${PYTHON}" "${RUNNER}" \
  --root "${ROOT}" \
  --runtime "${RUNTIME}" \
  --python "${PYTHON}" \
  --campaign "${CAMPAIGN}" \
  --molecular "${MOLECULAR}" \
  --folds "${FOLDS}" \
  --fit-views "${FIT_VIEWS}" \
  --preflight "${PREFLIGHT}" \
  --diagnostic "${DIAGNOSTIC}" \
  --surface "${SURFACE}" \
  --fitter "${FITTER}" \
  --aggregator "${AGGREGATOR}" \
  --source "${RUNNER}" \
  --source "${AGGREGATOR}" \
  --source "${FITTER}" \
  --source "${SURFACE}" \
  --source "${ROOT}/src/masld_bench/artifacts.py" \
  --source "${ROOT}/src/masld_bench/hashing.py" \
  --source "${ROOT}/slurm/run_gse267145_histology_production_bundle_cpu.sbatch" \
  --fold-workers 5 \
  --blas-threads 3
