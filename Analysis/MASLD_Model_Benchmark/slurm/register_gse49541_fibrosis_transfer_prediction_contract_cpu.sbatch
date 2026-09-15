#!/usr/bin/env bash
#SBATCH --job-name=model-data-870
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
SOURCE_FIT=${ROOT}/executions/model-fit-869-21130053/source_fit
SOURCE_FIT_ARTIFACTS_SHA256=28bf82a9a7aace04ecf6ad550fa75cec7b39fc4d2c40a35781e5c7b02d665d10
ACTIVATION=${ROOT}/executions/model-data-860-21129432/activation
SUMMARY=${ROOT}/executions/model-data-086-21110619-summarization/summary
CONTRACT=${ROOT}/config/evaluation/gse49541_fibrosis_transfer_prediction_v1.json
TARGET=${ROOT}/executions/model-data-870-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" || -e "${CONTRACT}" ]]; then
  printf 'Refusing to overwrite the prediction contract or its registration.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${SOURCE_FIT_ARTIFACTS_SHA256}" "${SOURCE_FIT}/ARTIFACTS.json" \
  | sha256sum --check --strict
bash -n slurm/register_gse49541_fibrosis_transfer_prediction_contract_cpu.sbatch
"${PYTHON}" - "${SOURCE_FIT}" "${ACTIVATION}" "${SUMMARY}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
print("input frozen-tree verification passed")
PY
"${PYTHON}" scripts/register_gse49541_fibrosis_transfer_prediction_contract.py \
  --benchmark-root "${ROOT}" \
  --source-fit "${SOURCE_FIT}" \
  --source-fit-contract "${ROOT}/config/evaluation/gse267145_fibrosis_transfer_source_fit_v4.json" \
  --activation "${ACTIVATION}" \
  --external-expression "${SUMMARY}" \
  --output "${CONTRACT}" \
  > "${STAGE}/register_stdout.json" 2> "${STAGE}/register_stderr.txt"
cp -p "${CONTRACT}" "${STAGE}/gse49541_fibrosis_transfer_prediction_v1.json"
sha256sum "${CONTRACT}" > "${STAGE}/contract.sha256"
sha256sum \
  scripts/register_gse49541_fibrosis_transfer_prediction_contract.py \
  slurm/register_gse49541_fibrosis_transfer_prediction_contract_cpu.sbatch \
  config/evaluation/gse267145_fibrosis_transfer_source_fit_v4.json \
  "${SOURCE_FIT}/ARTIFACTS.json" \
  "${ACTIVATION}/ARTIFACTS.json" \
  "${SUMMARY}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
contract = json.loads(
    (root / "gse49541_fibrosis_transfer_prediction_v1.json").read_text(encoding="utf-8")
)
if (
    contract.get("status") != "registered_prediction_only_pending"
    or contract.get("external_labels_read") is not False
    or contract.get("external_metrics_allowed") is not False
    or not contract.get("fitted_model_ids")
):
    raise SystemExit("registered prediction contract differs")
freeze_tree(root, {
    "artifact_class": "gse49541_fibrosis_transfer_prediction_contract_registration",
    "selected_source_model_id": contract["selected_source_model_id"],
    "fitted_models": len(contract["fitted_model_ids"]),
    "unfittable_models": len(contract["unfittable_model_ids"]),
    "external_labels_read": False,
    "status": "passed",
})
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
cat "${TARGET}/contract.sha256"
