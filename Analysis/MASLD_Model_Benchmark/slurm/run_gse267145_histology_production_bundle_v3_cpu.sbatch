#!/usr/bin/env bash
#SBATCH --job-name=model-training-073
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=3
export MKL_NUM_THREADS=3
export OPENBLAS_NUM_THREADS=3
: "${BUNDLE_VALIDATION_ROOT:?submit with the frozen bundle-validation root export}"
: "${BUNDLE_VALIDATION_SHA256:?submit with the frozen bundle-validation SHA-256 export}"
export BUNDLE_VALIDATION_ROOT
export BUNDLE_VALIDATION_SHA256

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
CAMPAIGN=${ROOT}/executions/model-training-069-production-v3
MOLECULAR=${ROOT}/executions/model-data-064-21079902/fixture/molecular
FOLDS=${ROOT}/executions/model-data-064-21079902/fixture/folds
PREFLIGHT=${ROOT}/executions/model-training-067-21080776
DIAGNOSTIC=${ROOT}/executions/model-check-072-21088262
FIT_VIEWS=${PREFLIGHT}/fit_views
SURFACE=${ROOT}/config/evaluation/gse267145_histology_benchmark_surface.json
CONTRACT=${ROOT}/config/evaluation/gse267145_histology_production_v3.json
REUSE_MANIFEST=${ROOT}/config/evaluation/gse267145_histology_production_v3_reuse_manifest.json
FITTER=${ROOT}/scripts/fit_gse267145_histology_baselines_v3.py
BASE_FITTER=${ROOT}/scripts/fit_gse267145_histology_baselines.py
AGGREGATOR=${ROOT}/scripts/aggregate_gse267145_histology_production_predictions.py
RUNNER=${ROOT}/scripts/run_gse267145_histology_production_bundle_v3.py
BASE_RUNNER=${ROOT}/scripts/run_gse267145_histology_production_bundle.py

cd "${ROOT}"
PYTHONPATH="${ROOT}/src:${ROOT}" "${PYTHON}" "${RUNNER}" \
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
  --source "${BASE_RUNNER}" \
  --source "${AGGREGATOR}" \
  --source "${FITTER}" \
  --source "${BASE_FITTER}" \
  --source "${CONTRACT}" \
  --source "${REUSE_MANIFEST}" \
  --source "${SURFACE}" \
  --source "${ROOT}/src/masld_bench/artifacts.py" \
  --source "${ROOT}/src/masld_bench/hashing.py" \
  --source "${ROOT}/slurm/run_gse267145_histology_production_bundle_v3_cpu.sbatch" \
  --fold-workers 5 \
  --blas-threads 3
