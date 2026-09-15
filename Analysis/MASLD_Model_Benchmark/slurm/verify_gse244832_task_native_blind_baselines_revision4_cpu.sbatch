#!/usr/bin/env bash
#SBATCH --job-name=model-check-346
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
: "${SOURCE_JOB_ID:?SOURCE_JOB_ID must identify the completed production bundle}"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG=${ROOT}/config/campaigns/gse244832_task_native_blind_baselines_revision4_20260825.json
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SOURCE=${ROOT}/executions/model-data-345-${SOURCE_JOB_ID}
RUN=${ROOT}/executions/model-train-345-${SOURCE_JOB_ID}
TARGET=${ROOT}/executions/model-check-346-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite a baseline verification artifact.\n' >&2
  exit 1
fi
cd "${ROOT}"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
sha256sum -c <<'EOF'
4b82d55552f4feaf2f1902f3617092f3695f44898ba410733302600848e9761f  config/campaigns/gse244832_task_native_blind_baselines_revision4_20260825.json
ad68dadfbd84963a6273d3f3312fb2f2ab45e269eacbda9d48e87f2e61371dc9  src/masld_bench/gse244832_task_native_baselines.py
12e791300d53430c6189b5bb4e6486b014f28ab76ca38861ba5b09648ec51e35  scripts/run_gse244832_task_native_blind_baselines.py
cd5eead3b56907d648cb4b20fcffbd8f6611ef73fbd9d180508668c2a342d8f4  scripts/verify_gse244832_task_native_blind_baselines.py
29f913b8736012462dddc92d72304354ee9c4c92b1ec0941a601ef097e09bbee  scripts/commit_gse244832_atac_transport_predictions.py
1cbcde8af1a8cd530ff04d93ac55822916cf958d9b3b7f990ca5b7237580c806  tests/unit/test_gse244832_task_native_baselines.py
1ce3466a6f42a8d6de2c7cc7915e4e0ef089a46878da3aa7de44b3e8205a1ede  tests/unit/test_run_gse244832_task_native_blind_baselines.py
85a71ccaa5fb7d67ff3d92319918f4943e494a3f521ab9c51a413ced4ca63546  tests/unit/test_verify_gse244832_task_native_blind_baselines.py
194a0a9d96ce036ec1e67409dfd59623c079b706ca3547accf67d9f9278c0a6f  executions/model-check-337-21101061/ARTIFACTS.json
EOF
PYTHONPATH=src:. "${PYTHON}" -m unittest \
  tests.unit.test_gse244832_task_native_baselines \
  tests.unit.test_run_gse244832_task_native_blind_baselines \
  tests.unit.test_verify_gse244832_task_native_blind_baselines \
  tests.unit.test_commit_gse244832_atac_transport_predictions
SOURCE_SHA=$(sha256sum "${SOURCE}/ARTIFACTS.json" | cut -d' ' -f1)
RUN_SHA=$(sha256sum "${RUN}/ARTIFACTS.json" | cut -d' ' -f1)
PYTHONPATH=src:. "${PYTHON}" scripts/verify_gse244832_task_native_blind_baselines.py \
  --root "${ROOT}" --config "${CONFIG}" --source-input "${SOURCE}" \
  --source-artifacts-sha256 "${SOURCE_SHA}" --run-input "${RUN}" \
  --run-artifacts-sha256 "${RUN_SHA}" --output "${STAGE}"
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
