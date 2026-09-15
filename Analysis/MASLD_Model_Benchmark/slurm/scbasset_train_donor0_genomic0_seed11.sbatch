#!/usr/bin/env bash
#SBATCH --job-name=masld-scbasset-train
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${OUTER_FOLD:=0}"
: "${SEED:=20260824}"
: "${INPUT_JOB_ID:=21064455}"
: "${INPUT_ARTIFACTS_SHA256:=02f76f4f7af62063a606943e2b1bcf11ca0c8d8a63f29879ed31bfc1fbc30dc0}"
: "${INPUT_ROOT_OVERRIDE:=}"
case "${OUTER_FOLD}" in 0|1|2|3|4) ;; *) printf 'Invalid OUTER_FOLD\n' >&2; exit 2 ;; esac
case "${SEED}" in 20260824|20260825|20260826) ;; *) printf 'Invalid screening SEED\n' >&2; exit 2 ;; esac
export PYTHONHASHSEED="${SEED}"
SPLIT_ID="donor${OUTER_FOLD}_genomic${OUTER_FOLD}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
if [[ -n "${INPUT_ROOT_OVERRIDE}" ]]; then
  case "${INPUT_ROOT_OVERRIDE}" in
    "${ROOT}/executions/scbasset-donor0-genomic0-inputs-21064455"|\
    "${ROOT}/executions/scbasset-donor1_genomic1-inputs-21067171"|\
    "${ROOT}/executions/scbasset-donor2_genomic2-inputs-21067172"|\
    "${ROOT}/executions/scbasset-donor3_genomic3-inputs-21067173"|\
    "${ROOT}/executions/scbasset-donor4_genomic4-inputs-21067174") ;;
    *) printf 'Invalid INPUT_ROOT_OVERRIDE\n' >&2; exit 2 ;;
  esac
  INPUT_ROOT="${INPUT_ROOT_OVERRIDE}"
else
  INPUT_ROOT="${ROOT}/executions/scbasset-${SPLIT_ID}-inputs-${INPUT_JOB_ID}"
fi
INPUTS="${INPUT_ROOT}/inputs"
if [[ "${INPUT_ARTIFACTS_SHA256}" == AUTO ]]; then
  if [[ ! -f "${INPUT_ROOT}/ARTIFACTS.json" ]]; then
    printf 'Dependency-complete input receipt is missing: %s\n' "${INPUT_ROOT}/ARTIFACTS.json" >&2
    exit 1
  fi
  read -r INPUT_ARTIFACTS_SHA256 _ < <(sha256sum "${INPUT_ROOT}/ARTIFACTS.json")
fi
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/scbasset-aed3a6f713091fd988196b297e9c06e092ff1d22-21064184"
SOURCE="${ACQUISITION}/source/scBasset-aed3a6f713091fd988196b297e9c06e092ff1d22"
SOURCE_ARTIFACTS_SHA256=9b0b418cea67bc788d7657913965c5e11ed62e36c35a8104c381c9007f6dffc5
RUNTIME_PROBE="${ROOT}/executions/scbasset-native-runtime-probe-21064213"
TRAINING_PROBE="${ROOT}/executions/scbasset-training-script-probe-21064318"
TARGET="${ROOT}/executions/scbasset-${SPLIT_ID}-seed${SEED}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" ]]; then
    if [[ -e "${failed}" ]]; then
      printf 'Cannot preserve failed scBasset training stage: %s\n' "${failed}" >&2
      return
    fi
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed scBasset training stage: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite scBasset training campaign: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

printf '%s  %s\n' "${INPUT_ARTIFACTS_SHA256}" "${INPUT_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" "${ACQUISITION}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' a46198d7e60a8f88dabdb201aed0f3738cc823909d3c9ba2dd63ceaa75668844 \
  "${RUNTIME_PROBE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 6e8cc84027ea4ce97d3f3b6dce0a3efe947c4b22be83060a88ac7412b8fd06c4 \
  "${TRAINING_PROBE}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${INPUT_ROOT}" "${BASE}" "${CORE}" \
  "${ACQUISITION}" "${RUNTIME_PROBE}" "${TRAINING_PROBE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("scBasset production training inputs verified")
PY

module purge
module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv,noheader \
  > "${STAGE}/validation/nvidia-smi.csv"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}/src:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python "${ROOT}/scripts/scbasset_train_fold.py" \
    --inputs "${INPUTS}" \
    --source "${SOURCE}" \
    --output "${STAGE}/model" \
    --input-artifacts-sha256 "${INPUT_ARTIFACTS_SHA256}" \
    --source-artifacts-sha256 "${SOURCE_ARTIFACTS_SHA256}" \
    --split-id "${SPLIT_ID}" \
    --donor-test-fold "${OUTER_FOLD}" \
    --donor-valid-fold "$(( (OUTER_FOLD + 1) % 5 ))" \
    --seed "${SEED}" \
    --batch-size 64 \
    --learning-rate 0.01 \
    --max-epochs 100 \
    --early-stopping-patience 10 \
    > "${STAGE}/validation/training.stdout.txt"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:ro" \
  --env "PYTHONPATH=${ROOT}/src:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/src/masld_bench/adapters/scbasset_native.py" \
  "${ROOT}/scripts/scbasset_train_fold.py" \
  "${ROOT}/slurm/scbasset_train_donor0_genomic0_seed11.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${INPUT_ARTIFACTS_SHA256}" \
  "${INPUT_ROOT}" "${SPLIT_ID}" "${SEED}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
input_root = Path(sys.argv[3]).resolve(strict=True)
split_id = sys.argv[4]
seed = int(sys.argv[5])
receipt = json.loads((root / "model" / "model_receipt.json").read_text())
input_summary = json.loads((input_root / "inputs" / "summary.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["split_id"] != split_id
    or receipt["donor_test_fold"] != int(split_id.split("_", 1)[0].removeprefix("donor"))
    or receipt["donor_valid_fold"] != (int(split_id.split("_", 1)[0].removeprefix("donor")) + 1) % 5
    or receipt["seed"] != seed
    or receipt["n_training_cells"] != input_summary["training_cells"]
    or receipt["public_tutorial_weights_used"]
    or receipt["held_donor_atac_used"]
    or receipt["genomic_test_atac_used"]
    or receipt["genomic_test_atac_available_to_training"]
):
    raise SystemExit("scBasset production model receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "scbasset_fold_model",
        "model_id": "scbasset",
        "dataset_id": "gse296875",
        "split_id": split_id,
        "seed": seed,
        "input_artifacts_sha256": sys.argv[2],
        "source_artifacts_sha256": "9b0b418cea67bc788d7657913965c5e11ed62e36c35a8104c381c9007f6dffc5",
        "runtime_probe_artifacts_sha256": "a46198d7e60a8f88dabdb201aed0f3738cc823909d3c9ba2dd63ceaa75668844",
        "training_initialization": "from_scratch",
        "public_tutorial_weights_used": False,
        "held_donor_atac_used": False,
        "genomic_test_atac_used": False,
        "genomic_test_atac_available_to_training": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scBasset %s seed%s model: %s\n' "${SPLIT_ID}" "${SEED}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
