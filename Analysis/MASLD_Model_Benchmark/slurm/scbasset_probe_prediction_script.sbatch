#!/usr/bin/env bash
#SBATCH --job-name=masld-scbasset-predict-probe
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/scbasset-aed3a6f713091fd988196b297e9c06e092ff1d22-21064184"
SOURCE="${ACQUISITION}/source/scBasset-aed3a6f713091fd988196b297e9c06e092ff1d22"
TRAINING_PROBE="${ROOT}/executions/scbasset-training-script-probe-21064249"
TARGET="${ROOT}/executions/scbasset-prediction-script-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite scBasset prediction probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/inputs" "${STAGE}/validation"
printf '%s  %s\n' 767c09734bd6c8c5021a37a6477f331f4d5d28c63019e8150097620788961945 \
  "${TRAINING_PROBE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 9b0b418cea67bc788d7657913965c5e11ed62e36c35a8104c381c9007f6dffc5 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict

cp --reflink=auto "${TRAINING_PROBE}/inputs/valid_seqs.h5" "${STAGE}/inputs/valid_seqs.h5"
cp --reflink=auto "${TRAINING_PROBE}/inputs/training_cells.tsv" "${STAGE}/inputs/training_cells.tsv"
module load python/3.11.5-GCCcore-13.2.0
python - "${STAGE}/inputs" <<'PY'
from pathlib import Path
import sys

root = Path(sys.argv[1])
(root / "valid_regions.tsv").write_text(
    "region_id\tblock_id\trole\tsequence_start\tsequence_end\n"
    "valid_1\tchr2\tvalid\t1000\t2344\n"
    "valid_2\tchr8\tvalid\t3000\t4344\n"
    "valid_3\tchr13\tvalid\t5000\t6344\n"
    "valid_4\tchr19\tvalid\t7000\t8344\n",
    encoding="utf-8",
)
(root / "held_valid_donors.tsv").write_text(
    "donor_id\touter_fold\tevaluation_role\n"
    "held_a\t1\tvalid\n"
    "held_b\t1\tvalid\n",
    encoding="utf-8",
)
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
  "${IMAGE}" python "${ROOT}/scripts/scbasset_predict_fold.py" \
    --inputs "${STAGE}/inputs" \
    --model "${TRAINING_PROBE}/model" \
    --source "${SOURCE}" \
    --output "${STAGE}/predictions" \
    --role valid \
    --run-id cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc \
    --join-namespace scbasset:fixture:valid \
    --input-artifacts-sha256 aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa \
    --model-artifacts-sha256 767c09734bd6c8c5021a37a6477f331f4d5d28c63019e8150097620788961945 \
    --source-artifacts-sha256 9b0b418cea67bc788d7657913965c5e11ed62e36c35a8104c381c9007f6dffc5 \
    --batch-size 2 \
    > "${STAGE}/validation/prediction.stdout.txt"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "predictions" / "prediction_receipt.json").read_text())
bundle = json.loads((root / "predictions" / "bundle" / "prediction_bundle.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["training_cells"] != 15
    or receipt["regions"] != 4
    or receipt["held_donors"] != 2
    or receipt["held_donor_atac_used"]
    or receipt["held_cell_embedding_available"]
    or bundle["n_predictions"] != 40
    or not bundle["metadata"]["profile_fixture_passed"]
):
    raise SystemExit("scBasset prediction probe differs")
freeze_tree(
    root,
    {
        "artifact_class": "scbasset_prediction_script_probe",
        "model_id": "scbasset",
        "device_class": "NVIDIA_L40S_sm89",
        "reverse_complement_tta": True,
        "profile_fixture_passed": True,
        "held_donor_atac_used": False,
        "held_cell_embedding_available": False,
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scBasset prediction script probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
