#!/usr/bin/env bash
#SBATCH --job-name=model-data-061
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
JOIN=${ROOT}/executions/gse267145-authoritative-join-21064930
MEASUREMENT=${ROOT}/executions/gse267145-rna-measurement-21065097
REFERENCE=${ROOT}/executions/gse267145-reference-crosswalk-21065154
QC_RIGHTS=${ROOT}/executions/gse267145-qc-rights-21065336
COORDINATE_ACTIVATION=${ROOT}/executions/gse267145-coordinate-activation-21070593
TASK_SPEC=${ROOT}/config/evaluation/gse267145_histology_state_task.toml
PROMOTION_GATE=${ROOT}/config/evaluation/gse267145_histology_state_promotion_gate.json
TARGET=${ROOT}/executions/model-data-061-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite histology-task activation artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}/src:${ROOT}"

python -m unittest -v tests.unit.test_activate_gse267145_histology_state \
  > "${STAGE}/unit_tests.txt" 2>&1
python - "${JOIN}" "${MEASUREMENT}" "${REFERENCE}" "${QC_RIGHTS}" "${COORDINATE_ACTIVATION}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
print("input frozen-tree verification passed")
PY
python scripts/activate_gse267145_histology_state.py \
  --join "${JOIN}" \
  --join-artifacts-sha256 e6c539c5fb29cd357dc126073779948c73c27fd219ff7a4d80b3d20a6cddb580 \
  --measurement "${MEASUREMENT}" \
  --measurement-artifacts-sha256 38531e463837cafd2c8e321411048c7b8b493876c3c964c0f50790cac8806f2f \
  --reference "${REFERENCE}" \
  --reference-artifacts-sha256 da491302ea738a92017978518239bc3884a7635d580e01e3bd28cd4f55b630c2 \
  --qc-rights "${QC_RIGHTS}" \
  --qc-rights-artifacts-sha256 ac7c92fd01b2b4240b380e4f6a58d98f1a6f853f8d0ad98d7167d7af34baaa73 \
  --coordinate-activation "${COORDINATE_ACTIVATION}" \
  --coordinate-activation-artifacts-sha256 e5042eb43dddd84cda81c5aeed5d6dbeb0378ab1c942a0183e2ef1b41cd18b20 \
  --task-spec "${TASK_SPEC}" \
  --task-spec-sha256 285376918c196341158e5184393a04946206062894d6596a571874aeb29bc3d8 \
  --promotion-gate "${PROMOTION_GATE}" \
  --promotion-gate-sha256 c8fd739e7d9cbaca9bf64729a8748dd0add42a9c76a5f7d4088362caf5add494 \
  --output "${STAGE}/activation" \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
python -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/activate_gse267145_histology_state.py \
  tests/unit/test_activate_gse267145_histology_state.py \
  slurm/activate_gse267145_histology_state.sbatch \
  "${TASK_SPEC}" \
  "${PROMOTION_GATE}" \
  "${JOIN}/ARTIFACTS.json" \
  "${MEASUREMENT}/ARTIFACTS.json" \
  "${REFERENCE}/ARTIFACTS.json" \
  "${QC_RIGHTS}/ARTIFACTS.json" \
  "${COORDINATE_ACTIVATION}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
activation = root / "activation"
contract = json.loads((activation / "activation_contract.json").read_text())
summary = json.loads((activation / "endpoint_summary.json").read_text())
if (
    contract.get("task_specific_model_fitting_allowed") is not True
    or contract.get("claim_contract", {}).get("champion_eligible") is not False
    or summary.get("participants") != 99
):
    raise SystemExit("histology activation contract differs")
freeze_tree(
    activation,
    {
        "artifact_class": "gse267145_histology_state_development_view",
        "cohort_family_id": "gse267145_znf469_human_liver",
        "task_id": "paired_bulk_histology_state",
        "participants": 99,
        "biological_unit": "participant",
        "pairing": "same_sample_different_aliquot",
        "outcome_visibility": "development_only",
        "task_specific_model_fitting_allowed": True,
        "champion_eligible": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
verify_frozen_tree(activation)
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_histology_state_activation_campaign",
        "task_id": "paired_bulk_histology_state",
        "participants": 99,
        "model_fitted": False,
        "metrics_calculated": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/activation/ARTIFACTS.json"
