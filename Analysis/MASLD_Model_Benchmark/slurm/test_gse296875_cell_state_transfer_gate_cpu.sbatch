#!/usr/bin/env bash
#SBATCH --job-name=model-check-335
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
TARGET=${ROOT}/executions/model-check-335-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
trap 'status=$?; if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi; exit ${status}' EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then exit 1; fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="${ROOT}:${ROOT}/src"
"${RUNTIME}/env/bin/python" -m unittest -v \
  tests.unit.test_gse296875_cell_state_transfer_gate \
  tests.unit.test_gse296875_transcriptformer_gpu_wrapper \
  > "${STAGE}/unit_tests.txt" 2>&1
bash -n slurm/run_gse296875_transcriptformer_query_embeddings.sbatch
"${RUNTIME}/env/bin/python" -m py_compile \
  scripts/admit_gse296875_cell_state_transfer_inputs.py \
  scripts/build_gse296875_transcriptformer_query_fixture.py \
  scripts/predict_gse296875_cell_state_baseline_outer_ensemble.py \
  scripts/predict_gse296875_transcriptformer_outer_ensemble.py \
  scripts/lock_gse296875_cell_state_transfer_predictions.py \
  scripts/evaluate_gse296875_cell_state_transfer.py
sha256sum scripts/*gse296875*cell_state* scripts/*gse296875*transcriptformer* \
  tests/unit/test_gse296875_cell_state_transfer_gate.py \
  tests/unit/test_gse296875_transcriptformer_gpu_wrapper.py \
  slurm/run_gse296875_transcriptformer_query_embeddings.sbatch \
  > "${STAGE}/source.sha256"
"${RUNTIME}/env/bin/python" - <<'PY' "${STAGE}"
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "gse296875_cell_state_transfer_gate_cpu_validation",
    "tests": 7,
    "gpu_jobs_run": False,
    "evaluator_labels_read": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
