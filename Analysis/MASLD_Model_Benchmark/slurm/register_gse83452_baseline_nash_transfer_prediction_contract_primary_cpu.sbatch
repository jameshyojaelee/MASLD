#!/usr/bin/env bash
#SBATCH --job-name=model-data-884
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
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
ARM_ID=nash_vs_nafl
SOURCE_FIT=${ROOT}/executions/model-fit-882-21130285-nash_vs_nafl/source_fit
SOURCE_FIT_ARTIFACTS_SHA256=57d4f7c2505a064ac13aa038e0db798cc8370a3ef0816589ab8adab04e71f1d7
FIT_CONTRACT=${ROOT}/config/evaluation/gse135251_nash_transfer_source_fit_nash_vs_nafl_v1.json
ACTIVATION=${ROOT}/executions/model-data-881-21130281-activation/activation
SUMMARY=${ROOT}/executions/model-data-088-21112137-summarization/summary
CONTRACT=${ROOT}/config/evaluation/gse83452_baseline_nash_transfer_prediction_nash_vs_nafl_v1.json
TARGET=${ROOT}/executions/model-data-884-${SLURM_JOB_ID}-nash_vs_nafl
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" || -e "${CONTRACT}" ]]; then
  printf 'Refusing to overwrite the GSE83452 prediction contract artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
bash -n slurm/register_gse83452_baseline_nash_transfer_prediction_contract_primary_cpu.sbatch
printf '%s  %s\n' "${SOURCE_FIT_ARTIFACTS_SHA256}" "${SOURCE_FIT}/ARTIFACTS.json" \
  | sha256sum --check --strict
"${PYTHON}" - "${SOURCE_FIT}" "${ACTIVATION}" "${SUMMARY}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
print("input frozen-tree verification passed")
PY
"${PYTHON}" scripts/register_gse83452_baseline_nash_transfer_prediction_contract.py \
  --benchmark-root "${ROOT}" \
  --source-fit "${SOURCE_FIT}" \
  --source-fit-contract "${FIT_CONTRACT}" \
  --activation "${ACTIVATION}" \
  --external-expression "${SUMMARY}" \
  --baseline-participants 152 \
  --output "${CONTRACT}" \
  > "${STAGE}/register_stdout.json" 2> "${STAGE}/register_stderr.txt"
cp -p "${CONTRACT}" "${STAGE}/prediction_contract.json"
sha256sum "${CONTRACT}" > "${STAGE}/contract.sha256"
sha256sum \
  scripts/register_gse83452_baseline_nash_transfer_prediction_contract.py \
  slurm/register_gse83452_baseline_nash_transfer_prediction_contract_primary_cpu.sbatch \
  "${FIT_CONTRACT}" \
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
contract = json.loads((root / "prediction_contract.json").read_text(encoding="utf-8"))
if (
    contract.get("status") != "registered_prediction_only_pending"
    or contract.get("external_labels_read") is not False
    or contract.get("external_metrics_allowed") is not False
    or contract.get("timepoint") != "baseline"
    or "age_sex_logistic" not in contract.get("unfittable_model_ids", [])
    or not contract.get("fitted_model_ids")
    or contract["selected_source_model_id"] in contract["unfittable_model_ids"]
):
    raise SystemExit("registered prediction contract differs")
freeze_tree(root, {
    "artifact_class": "gse83452_baseline_nash_transfer_prediction_contract_registration",
    "arm_id": contract["arm_id"],
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
