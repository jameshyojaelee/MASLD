#!/usr/bin/env bash
#SBATCH --job-name=masld-epibert-shape-audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
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
IMAGE="${BASE}/tensorflow-tensorflow-2.13.0-gpu.sif"
TARGET="${ROOT}/executions/epibert-checkpoint-shape-audit-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed EpiBERT shape audit: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite EpiBERT shape audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' "${BASE_ARTIFACTS_SHA256}" "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BASE}" "${RUNTIME}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
module purge
module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python "${ROOT}/scripts/epibert_checkpoint_shape_audit.py" \
    --pretrained-model1 "${RUNTIME}/checkpoints/pretrained_model1/ckpt-45" \
    --pretrained-model2 "${RUNTIME}/checkpoints/pretrained_model2/ckpt-45" \
    --fine-tuned-rampage "${RUNTIME}/checkpoints/fine_tuned_rampage/ckpt-16" \
    --output "${STAGE}/audit" \
    > "${STAGE}/audit.stdout.json"
sha256sum \
  "${ROOT}/scripts/epibert_checkpoint_shape_audit.py" \
  "${ROOT}/slurm/epibert_checkpoint_shape_audit.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "audit/receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["tensorflow"] != "2.13.0"
    or set(receipt["models"]) != {"pretrained_model1", "pretrained_model2", "fine_tuned_rampage"}
    or receipt["constructor_changed_by_audit"]
    or receipt["project_data_read"]
    or receipt["outcomes_read"]
):
    raise SystemExit("EpiBERT checkpoint shape audit differs")
freeze_tree(root, {
    "artifact_class": "epibert_checkpoint_shape_audit",
    "model_roster": sorted(receipt["models"]),
    "project_data_read": False,
    "outcomes_read": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'EpiBERT checkpoint shape audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
