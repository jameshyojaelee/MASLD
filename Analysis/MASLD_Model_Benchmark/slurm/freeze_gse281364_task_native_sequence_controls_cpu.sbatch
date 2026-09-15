#!/usr/bin/env bash
#SBATCH --job-name=model-check-227
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
CONVERSION="${ROOT}/executions/environments/dna-lm-torch2.6-conversion-21064945"
CONFIG="${ROOT}/config/gse281364_task_native_sequence_controls.toml"
TARGET="${ROOT}/executions/model-check-227-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed sequence-control fixture: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite sequence-control fixture: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

for binding in \
  "61f3ce0d8a92df554ba4bfdbaa674d4838a4c0e3762025f3f06cb5c1fcf5d00e executions/model-check-220-21088696" \
  "33c3ea9fe134157b9f8f29396cd9d0d565243a9f6e44db3f458cd970c486d32c executions/gse281364-mpralegnet-fixture-21069059" \
  "fdc1db620152757f74c07904487a62381fcd182a532669b93b0aa1ca5cdfc2ac executions/model-check-225-21091492"; do
  read -r expected relative <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}/ARTIFACTS.json" \
    | sha256sum --check --strict
done

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
for relative in (
    "executions/model-check-220-21088696",
    "executions/gse281364-mpralegnet-fixture-21069059",
    "executions/model-check-225-21091492",
):
    verify_frozen_tree(root / relative)
print("task-native sequence-control authorities verified")
PY
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  "${ROOT}/scripts/build_gse281364_task_native_sequence_fixture.py" \
  --project-root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/frozen" \
  > "${STAGE}/validation/build.stdout.json" \
  2> "${STAGE}/validation/build.stderr.txt"
module purge

module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${CONVERSION}/site-packages:${TORCH}/site-packages:${BASE}/site-packages:${ROOT}/src:${ROOT}"
python -m unittest -v tests.unit.test_gse281364_task_native_sequence_controls \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
python "${ROOT}/scripts/validate_gse281364_task_native_sequence_controls.py" \
  --fixture "${STAGE}/frozen" \
  --architecture-source "${ROOT}/scripts/gse281364_task_native_sequence_controls.py" \
  --output "${STAGE}/frozen/validation" \
  > "${STAGE}/validation/architecture.stdout.json" \
  2> "${STAGE}/validation/architecture.stderr.txt"
python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${ROOT}/config/gse281364_task_native_sequence_controls.toml" \
  "${ROOT}/scripts/gse281364_task_native_sequence_controls.py" \
  "${ROOT}/scripts/build_gse281364_task_native_sequence_fixture.py" \
  "${ROOT}/scripts/validate_gse281364_task_native_sequence_controls.py" \
  "${ROOT}/tests/unit/test_gse281364_task_native_sequence_controls.py" \
  "${ROOT}/slurm/freeze_gse281364_task_native_sequence_controls_cpu.sbatch" \
  > "${STAGE}/source.sha256"
module purge

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
fixture = json.loads((root / "frozen/fixture/receipt.json").read_text())
validation = json.loads((root / "frozen/validation/architecture_contract.json").read_text())
task = json.loads((root / "frozen/contract/task_spec.json").read_text())
license_value = json.loads(
    (root / "frozen/contract/code_license_disposition.json").read_text()
)
if (
    fixture["status"] != "pass_outcome_blind_task_fixture_freeze"
    or validation["status"] != "pass_outcome_blind_architecture_fixture"
    or not validation["five_genuine_seed_states_per_model"]
    or not validation["shared_ref_alt_encoder"]
    or validation["local_atac_code_reused"]
    or validation["local_atac_weights_read"]
    or fixture["outcomes_read"]
    or validation["outcomes_read"]
    or fixture["metrics_calculated"]
    or validation["metrics_calculated"]
    or task["training_executed"]
    or task["predictions_generated"]
    or license_value["source_redistribution"]
    != "blocked_pending_repository_license"
):
    raise SystemExit("task-native sequence-control freeze differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse281364_task_native_sequence_control_task_and_fixture",
        "dataset_id": "gse281364",
        "models": ["sequence_cnn_control", "sequence_transformer_control"],
        "elements": 1033,
        "long_range_blocks": 239,
        "genuine_training_seeds": [1103, 2909, 4721, 6673, 8111],
        "outcomes_read": False,
        "metrics_calculated": False,
        "model_training_executed": False,
        "predictions_generated": False,
        "local_atac_weights_read": False,
        "release_eligible": False,
        "status": "passed",
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Task-native sequence-control freeze: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
