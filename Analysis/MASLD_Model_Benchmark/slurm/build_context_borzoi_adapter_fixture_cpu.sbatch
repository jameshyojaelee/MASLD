#!/usr/bin/env bash
#SBATCH --job-name=model-check-982
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONTRACT="${ROOT}/config/artifacts/models/context_borzoi/rna_context_adapter_fixture_20260825.json"
TARGET="${ROOT}/executions/model-check-982-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite adapter fixture: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_context_borzoi_adapter_fixture \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/build_context_borzoi_adapter_fixture.py \
  --project-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/fixture" \
  > "${STAGE}/validation/build.stdout.json"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  config/artifacts/models/context_borzoi/rna_context_adapter_fixture_20260825.json \
  config/artifacts/models/context_borzoi/component_interface_readiness_20260825.json \
  config/campaigns/v1_conditional_model.toml \
  executions/model-check-981-21099281/ARTIFACTS.json \
  src/masld_bench/context_preconditions.py \
  src/masld_bench/adapters/context_borzoi_fixture.py \
  scripts/build_context_borzoi_adapter_fixture.py \
  tests/unit/test_context_borzoi_adapter_fixture.py \
  slurm/build_context_borzoi_adapter_fixture_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.adapters.context_borzoi_fixture import verify_nonpickle_export

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "fixture/receipt.json").read_text(encoding="utf-8"))
stderr = (root / "validation/tests.stderr.txt").read_text(encoding="utf-8")
verify_nonpickle_export(root / "fixture/fixture_export")
if (
    receipt.get("status") != "pass_synthetic_fixture_architecture_blocked"
    or receipt.get("component_authority_verified") is not True
    or receipt.get("complementarity_trigger_passed") is not False
    or receipt.get("component_selection_locked") is not False
    or receipt.get("parent_component_selected") is not False
    or receipt.get("architecture_built") is not False
    or receipt.get("checkpoint_opened") is not False
    or receipt.get("model_weights_present") is not False
    or receipt.get("training_performed") is not False
    or receipt.get("biological_prediction_performed") is not False
    or receipt.get("biological_data_read") is not False
    or receipt.get("development_outcomes_read") is not False
    or receipt.get("sealed_features_or_outcomes_read") is not False
    or receipt.get("observed_atac_read") is not False
    or receipt.get("observed_multiome_task_modified") is not False
    or receipt.get("gpu_executed") is not False
    or receipt.get("sequence_only_exact_identity") is not True
    or receipt.get("missing_context_exact_sequence_only_identity") is not True
    or receipt.get("trans_only_sequence_invariant") is not True
    or receipt.get("permuted_context_deranged_within_stratum") is not True
    or receipt.get("pickle_present") is not False
    or receipt.get("tensor_state_present") is not False
    or receipt.get("conditional_model_authorized") is not False
    or receipt.get("gpu_job_submission_authorized") is not False
    or "Ran 9 tests" not in stderr
    or not stderr.rstrip().endswith("OK")
):
    raise SystemExit("context-Borzoi adapter fixture receipt differs")

freeze_tree(
    root,
    {
        "artifact_class": "context_borzoi_rna_adapter_synthetic_fixture",
        "model_id": "context_borzoi",
        "artifact_role": "biology_free_cpu_interface_fixture_only",
        "tests": 9,
        "component_authority_verified": True,
        "complementarity_trigger_passed": False,
        "component_selection_locked": False,
        "parent_component_selected": False,
        "architecture_built": False,
        "checkpoint_opened": False,
        "model_weights_present": False,
        "training_performed": False,
        "biological_prediction_performed": False,
        "biological_data_read": False,
        "development_outcomes_read": False,
        "sealed_features_or_outcomes_read": False,
        "observed_atac_read": False,
        "observed_multiome_task_modified": False,
        "gpu_executed": False,
        "export_format": "canonical_json_utf8",
        "pickle_present": False,
        "conditional_model_authorized": False,
        "status": "passed",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Context-Borzoi RNA adapter fixture: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
