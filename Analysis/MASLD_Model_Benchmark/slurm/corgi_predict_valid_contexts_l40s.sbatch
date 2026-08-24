#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-valid-predict
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${OUTER_FOLD:?OUTER_FOLD is required}"
: "${MAPPER_JOB_ID:?MAPPER_JOB_ID is required}"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/corgi-nvidia-pytorch-24.01-py3-recovery-21064690"
IMAGE="${RUNTIME}/nvidia-pytorch-24.01-py3.sif"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/corgi-d84c000-zenodo18630048-21064463"
SOURCE="${ACQUISITION}/source/corgi-d84c00082da3d04eb05169e290e53f1c7da2b9c5"
CHECKPOINT="${ACQUISITION}/downloads/corgi_model.pt"
OVERLAY="${ROOT}/executions/environments/corgi-rotary-emb-0.1-flashattn2.0.4-build-r1-21064845"
RUNTIME_PROBE="${ROOT}/executions/corgi-native-runtime-l40s-r8-21065022"
MAPPER="${ROOT}/executions/corgi-gse296875-mapper-fold${OUTER_FOLD}-${MAPPER_JOB_ID}"
FIXTURE="${ROOT}/executions/borzoi-enformer-sequence-fixture-21064429"
FASTA="${FIXTURE}/fixture/borzoi_ensemble.reference.fa.gz"
TARGET="${ROOT}/executions/corgi-gse296875-valid-prediction-fold${OUTER_FOLD}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi prediction artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in \
  "01de778c17c67d92de9a2e242350e17443dbfe76b7c5a5f8bc447eae5408cd3c ${RUNTIME}/ARTIFACTS.json" \
  "73aad9ab9ad3e080a2ccaa0cee1a1d727d126a454f0e69f74f343a900b168c0f ${ACQUISITION}/ARTIFACTS.json" \
  "c60c4ebf0b32df7d158a3a1d781f5d45fdf4187ab0fba9f4bea35605a92e4c3a ${OVERLAY}/ARTIFACTS.json" \
  "09d8811b612cba4f75361feeba29bd0676efde83f0e4785a37990c2d2f52a0a2 ${RUNTIME_PROBE}/ARTIFACTS.json" \
  "d8a44d36c5b2b9997971e358e3553c10ce24c4badd41b28c91cfd48986123f0e ${FIXTURE}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${MAPPER}" "${RUNTIME_PROBE}" "${FIXTURE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY
install -d -m 0750 "${STAGE}"
module purge
module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total --format=csv,noheader > "${STAGE}/nvidia-smi.csv"
singularity exec --nv --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env PYTHONHASHSEED=20260824 \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  --env "PYTHONPATH=${ROOT}:${OVERLAY}/site-packages" \
  "${IMAGE}" python "${ROOT}/tests/unit/test_corgi_predict_valid_contexts.py"
singularity exec --nv --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env PYTHONHASHSEED=20260824 \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  --env "PYTHONPATH=${ROOT}:${OVERLAY}/site-packages" \
  "${IMAGE}" python "${ROOT}/scripts/corgi_predict_valid_contexts.py" \
    --source "${SOURCE}" --checkpoint "${CHECKPOINT}" --mapper "${MAPPER}" \
    --fixture-fasta "${FASTA}" --batch-size 4 --output "${STAGE}/prediction"
singularity exec --nv --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" --env PYTHONNOUSERSITE=1 "${IMAGE}" python -m pip freeze > "${STAGE}/pip-freeze.txt"
sha256sum scripts/corgi_predict_valid_contexts.py scripts/corgi_probe_native_runtime.py tests/unit/test_corgi_predict_valid_contexts.py slurm/corgi_predict_valid_contexts_l40s.sbatch > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${OUTER_FOLD}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root, fold = Path(sys.argv[1]).resolve(strict=True), int(sys.argv[2])
receipt = json.loads((root / "prediction/receipt.json").read_text())
if receipt["status"] != "pass_outcome_free_prediction" or receipt["held_ATAC_or_other_outcomes_used"]:
    raise SystemExit("Corgi prediction receipt differs")
freeze_tree(root, {"artifact_class": "Corgi_valid_role_context_predictions", "model_id": "corgi_regular", "dataset_id": "gse296875", "outer_fold": fold, "records": receipt["records"], "outcomes_read": False, "open_champion_eligible": False, "status": "passed"})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi valid context predictions: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
