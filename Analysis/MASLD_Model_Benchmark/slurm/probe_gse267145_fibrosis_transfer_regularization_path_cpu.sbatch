#!/usr/bin/env bash
#SBATCH --job-name=model-check-867
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
CONTRACT=${ROOT}/config/evaluation/gse267145_fibrosis_transfer_source_fit_v2.json
TARGET=${ROOT}/executions/model-check-867-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite regularization-probe artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/probe_gse267145_fibrosis_transfer_regularization_path_cpu.sbatch
"${PYTHON}" scripts/probe_gse267145_fibrosis_transfer_regularization_path.py \
  --benchmark-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/regularization_path.json" \
  > "${STAGE}/probe_stdout.json" 2> "${STAGE}/probe_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/probe_gse267145_fibrosis_transfer_regularization_path.py \
  slurm/probe_gse267145_fibrosis_transfer_regularization_path_cpu.sbatch \
  "${CONTRACT}" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(root, {
    "artifact_class": "gse267145_fibrosis_transfer_regularization_probe",
    "external_expression_values_read": False,
    "external_labels_read": False,
    "metrics_calculated": False,
    "status": "passed",
})
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
cat "${TARGET}/probe_stdout.json"
