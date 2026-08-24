#!/usr/bin/env bash
#SBATCH --job-name=masld-epibert-released-forward
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/epibert-tensorflow-2.13.0-21064371"
BASE_ARTIFACTS_SHA256=10a8691e9468dfa44a25e55770b477e77d2ebed12839d49203f1a386596af27f
RUNTIME="${ROOT}/executions/environments/epibert-executable-tf2.13-21064522"
RUNTIME_ARTIFACTS_SHA256=83e42afe350081391204703e9178ac30b52b7dfc15b863ba1af2e10611342aa4
FIXTURE="${ROOT}/executions/epibert-epcotv2-no-outcome-fixture-21064487"
FIXTURE_ARTIFACTS_SHA256=4ac1e23e884dab0318c2c2fbf18e5a8ea244c5db8d0bf40ef95f3479b4f25607
AUDIT="${ROOT}/executions/epibert-checkpoint-shape-audit-21065291"
AUDIT_ARTIFACTS_SHA256=937e455dcf7bd96081e8bb148925c0fe9b3b286a6a97529c6d48fa31287a6c12
IMAGE="${BASE}/tensorflow-tensorflow-2.13.0-gpu.sif"
SOURCE="${RUNTIME}/source/EpiBERT-dae61b434f885991c29e2fde3fbf7f70d0ff9b80"
SITE_PACKAGES="${RUNTIME}/site-packages"
INPUTS="${FIXTURE}/fixtures/epibert/inputs.npz"
TARGET="${ROOT}/executions/epibert-released-checkpoint-forward-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" ]]; then
    if [[ -e "${failed}" ]]; then
      printf 'Cannot preserve failed EpiBERT stage because target exists: %s\n' "${failed}" >&2
      return
    fi
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed EpiBERT stage: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite EpiBERT forward artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/models" "${STAGE}/validation"

printf '%s  %s\n' "${BASE_ARTIFACTS_SHA256}" "${BASE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${FIXTURE_ARTIFACTS_SHA256}" "${FIXTURE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${AUDIT_ARTIFACTS_SHA256}" "${AUDIT}/ARTIFACTS.json" \
  | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BASE}" "${RUNTIME}" "${FIXTURE}" "${AUDIT}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY
module purge
module load singularity/4.2.2

singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${SITE_PACKAGES}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m unittest -v \
    tests.unit.epibert_test_released_checkpoint_probe \
    > "${STAGE}/validation/tests.stdout.txt" \
    2> "${STAGE}/validation/tests.stderr.txt"

run_probe() {
  local name="$1"
  local kind="$2"
  local checkpoint="$3"
  local seed="$4"
  singularity exec --nv --cleanenv \
    --bind "${ROOT}:${ROOT}:rw" \
    --env "PYTHONPATH=${SOURCE}:${SITE_PACKAGES}" \
    --env PYTHONNOUSERSITE=1 \
    --env PYTHONDONTWRITEBYTECODE=1 \
    "${IMAGE}" python "${ROOT}/scripts/epibert_probe_released_checkpoint.py" \
      --model-kind "${kind}" \
      --checkpoint-prefix "${checkpoint}" \
      --fixture "${INPUTS}" \
      --output "${STAGE}/models/${name}" \
      --seed "${seed}" \
      > "${STAGE}/validation/${name}.stdout.json"
}

run_probe \
  pretrained_model1 \
  pretrained_atac \
  "${RUNTIME}/checkpoints/pretrained_model1/ckpt-45" \
  21
run_probe \
  pretrained_model2 \
  pretrained_atac \
  "${RUNTIME}/checkpoints/pretrained_model2/ckpt-45" \
  21
run_probe \
  fine_tuned_rampage \
  fine_tuned_rampage \
  "${RUNTIME}/checkpoints/fine_tuned_rampage/ckpt-16" \
  19

singularity exec --cleanenv \
  --env "PYTHONPATH=${SITE_PACKAGES}" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/epibert_probe_released_checkpoint.py" \
  "${ROOT}/tests/unit/epibert_test_released_checkpoint_probe.py" \
  "${ROOT}/slurm/epibert_probe_released_checkpoints.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${STAGE}" \
  "${BASE_ARTIFACTS_SHA256}" \
  "${RUNTIME_ARTIFACTS_SHA256}" \
  "${FIXTURE_ARTIFACTS_SHA256}" \
  "${AUDIT_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
names = ("pretrained_model1", "pretrained_model2", "fine_tuned_rampage")
receipts = {
    name: json.loads((root / f"models/{name}/receipt.json").read_text())
    for name in names
}
if set(receipts) != set(names):
    raise SystemExit("released EpiBERT roster differs")
for name, receipt in receipts.items():
    if (
        receipt["status"] != "pass"
        or receipt["tensorflow"] != "2.13.0"
        or receipt["mixed_precision_policy"] != "mixed_bfloat16"
        or "L40S" not in receipt["gpu"]
        or not receipt["restore_existing_objects_matched"]
        or not receipt["constructor_heads_bound_to_checkpoint_shape_audit"]
        or receipt["released_checkpoint_inferred_num_heads"] != 8
        or receipt["published_command_num_heads"]
        != (4 if name.startswith("pretrained_") else 8)
        or not receipt["repeat_inference_bit_identical"]
        or receipt["backward_executed"]
        or receipt["checkpoint_resume_executed"]
        or receipt["project_data_read"]
        or receipt["outcomes_read"]
        or receipt["labels_read"]
        or not receipt["observed_atac_required_at_query"]
        or not receipt["scored_atac_span_removed_from_input"]
        or receipt["rna_conditioned_atac_eligible"]
        or receipt["sealed_inference_eligible"]
    ):
        raise SystemExit(f"released EpiBERT receipt differs: {name}")
freeze_tree(
    root,
    {
        "artifact_class": "epibert_released_checkpoint_forward_probe",
        "model_roster": list(names),
        "base_runtime_artifacts_sha256": sys.argv[2],
        "executable_runtime_artifacts_sha256": sys.argv[3],
        "no_outcome_fixture_artifacts_sha256": sys.argv[4],
        "checkpoint_shape_audit_artifacts_sha256": sys.argv[5],
        "released_checkpoints_restored": True,
        "repeat_inference_bit_identical": True,
        "backward_executed": False,
        "checkpoint_resume_executed": False,
        "project_data_read": False,
        "outcomes_read": False,
        "rna_conditioned_atac_models": [],
        "sealed_inference_eligible_models": [],
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'EpiBERT released-checkpoint forward artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
