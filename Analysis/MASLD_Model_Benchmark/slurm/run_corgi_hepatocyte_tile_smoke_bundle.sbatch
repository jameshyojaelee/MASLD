#!/usr/bin/env bash
#SBATCH --job-name=model-work-101
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/corgi-nvidia-pytorch-24.01-py3-recovery-21064690"
IMAGE="${RUNTIME}/nvidia-pytorch-24.01-py3.sif"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/corgi-d84c000-zenodo18630048-21064463"
SOURCE="${ACQUISITION}/source/corgi-d84c00082da3d04eb05169e290e53f1c7da2b9c5"
CHECKPOINT="${ACQUISITION}/downloads/corgi_model.pt"
OVERLAY="${ROOT}/executions/environments/corgi-rotary-emb-0.1-flashattn2.0.4-build-r1-21064845"
RUNTIME_PROBE="${ROOT}/executions/corgi-native-runtime-l40s-r8-21065022"
TILE_CONTRACT="${ROOT}/executions/corgi-outcome-aligned-tiles-21068311"
TILE_ROSTER="${ROOT}/executions/corgi-tile-smoke-subset-21068454"
PREDICTOR_TEST="${ROOT}/executions/corgi-outcome-aligned-predictor-test-21068997"
TARGET="${ROOT}/executions/corgi-hepatocyte-tile-smoke-bundle-v3-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi bundle artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in \
  "01de778c17c67d92de9a2e242350e17443dbfe76b7c5a5f8bc447eae5408cd3c ${RUNTIME}/ARTIFACTS.json" \
  "73aad9ab9ad3e080a2ccaa0cee1a1d727d126a454f0e69f74f343a900b168c0f ${ACQUISITION}/ARTIFACTS.json" \
  "c60c4ebf0b32df7d158a3a1d781f5d45fdf4187ab0fba9f4bea35605a92e4c3a ${OVERLAY}/ARTIFACTS.json" \
  "09d8811b612cba4f75361feeba29bd0676efde83f0e4785a37990c2d2f52a0a2 ${RUNTIME_PROBE}/ARTIFACTS.json" \
  "1cdf94aa4c37ead139125ffc3d8fee358605a2394a206b1f3e79190ecd25645c ${TILE_CONTRACT}/ARTIFACTS.json" \
  "2617f3a85a9ced8b0f7cb9efe3819b1e57177629221a2013e39fc7c96e187f5b ${TILE_ROSTER}/ARTIFACTS.json" \
  "2639ccd153b854f69a06abb36f1456f1a1145677ddf18aaed5729509a4cc2f41 ${PREDICTOR_TEST}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
for pair in \
  "0c2f6577c9981039a39409d5fb3593973e089fd9d2d3d665f72947c55d1f7c7a ${ROOT}/executions/corgi-gse296875-mapper-fold0-21066334/ARTIFACTS.json" \
  "f964222d36d11755e44e77c07a1769f0ffebd3771dbd877909c3c2316911cc5a ${ROOT}/executions/corgi-gse296875-mapper-fold1-21066335/ARTIFACTS.json" \
  "ceae9d4a53ed6a8bcfaf08aeffe6a396aad12ee1839ebe2da818a75332709507 ${ROOT}/executions/corgi-gse296875-mapper-fold2-21066337/ARTIFACTS.json" \
  "18ddf520a41db2b968a766ced416c85c02bef355becef1c12f5c9b808052b367 ${ROOT}/executions/corgi-gse296875-mapper-fold3-21066336/ARTIFACTS.json" \
  "a2ef2ae490c3cfe429c0c29755d33285085adbe8bf36418e938345c31c3e6053 ${ROOT}/executions/corgi-gse296875-mapper-fold4-21066338/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
for pair in \
  "2a68bc359f246f4fecc6e835fc6fd798861b4c631fd9d665fccd85d200b6e31c ${ROOT}/scripts/corgi_predict_outcome_aligned_tiles.py" \
  "e351024065408c644b0476fffa41a002e6fd57f7b768972d020b3275fcfcf9f5 ${ROOT}/scripts/corgi_predict_valid_contexts.py" \
  "66c76174c9b5c0c9b0c55706176d2b6e74833476183c8682f1523be8f72fbc94 ${ROOT}/tests/unit/test_corgi_predict_outcome_aligned_tiles.py" \
  "f21d31c0c9888e3a80f730ff57ba4ccaf6b57cad887d542d16a442a241f34be6 ${ROOT}/scripts/corgi_probe_native_runtime.py"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${TILE_CONTRACT}" "${TILE_ROSTER}" "${PREDICTOR_TEST}" \
  "${ROOT}/executions/corgi-gse296875-mapper-fold0-21066334" \
  "${ROOT}/executions/corgi-gse296875-mapper-fold1-21066335" \
  "${ROOT}/executions/corgi-gse296875-mapper-fold2-21066337" \
  "${ROOT}/executions/corgi-gse296875-mapper-fold3-21066336" \
  "${ROOT}/executions/corgi-gse296875-mapper-fold4-21066338" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
for fold, value in enumerate(sys.argv[4:]):
    receipt = json.loads((Path(value) / "receipt.json").read_text())
    if (
        receipt["schema_version"] != "masld-bench-corgi-masked-context-mapper-v1"
        or receipt["status"] != "pass_outcome_free_mapper"
        or receipt["outer_fold"] != fold
        or receipt["valid_fold"] != (fold + 1) % 5
        or receipt["held_ATAC_or_other_outcomes_used"]
        or receipt["test_or_sealed_outcomes_read"]
    ):
        raise SystemExit(f"Corgi mapper fold {fold} receipt differs")
PY
if [[ "${PREFLIGHT_ONLY:-0}" == 1 ]]; then
  printf 'Corgi hepatocyte tile bundle preflight passed\n'
  exit 0
fi
install -d -m 0750 "${STAGE}"
module purge
module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total --format=csv,noheader > "${STAGE}/nvidia-smi.csv"
singularity exec --nv --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" --bind "${STAGE}:${STAGE}:rw" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  --env PYTHONHASHSEED=20260824 --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  --env "PYTHONPATH=${ROOT}:${OVERLAY}/site-packages" \
  "${IMAGE}" python -m unittest tests.unit.test_corgi_predict_outcome_aligned_tiles
MAPPER_JOBS=(21066334 21066335 21066337 21066336 21066338)
for fold in 0 1 2 3 4; do
  mapper="${ROOT}/executions/corgi-gse296875-mapper-fold${fold}-${MAPPER_JOBS[${fold}]}"
  singularity exec --nv --cleanenv --containall --net --network none \
    --bind "${ROOT}:${ROOT}:ro" --bind "${STAGE}:${STAGE}:rw" \
    --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
    --env PYTHONHASHSEED=20260824 --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
    --env "PYTHONPATH=${ROOT}:${OVERLAY}/site-packages" \
    "${IMAGE}" python "${ROOT}/scripts/corgi_predict_outcome_aligned_tiles.py" \
      --source "${SOURCE}" --checkpoint "${CHECKPOINT}" --mapper "${mapper}" \
      --tiles "${TILE_CONTRACT}/tiles" --tile-roster "${TILE_ROSTER}/subset" \
      --mapper-outer-fold "${fold}" --lineage hepatocyte --batch-size 4 \
      --output "${STAGE}/fold${fold}"
done
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/corgi_predict_outcome_aligned_tiles.py" \
  "${ROOT}/scripts/corgi_predict_valid_contexts.py" \
  "${ROOT}/scripts/corgi_probe_native_runtime.py" \
  "${ROOT}/tests/unit/test_corgi_predict_outcome_aligned_tiles.py" \
  "${ROOT}/slurm/run_corgi_hepatocyte_tile_smoke_bundle.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
for fold in range(5):
    receipt = json.loads((root / f"fold{fold}" / "receipt.json").read_text())
    valid_fold = (fold + 1) % 5
    if (
        receipt["schema_version"] != "masld-bench-corgi-outcome-aligned-prediction-v3"
        or receipt["status"] != "pass_outcome_free_prediction"
        or receipt["split_id"] != f"donor{fold}_genomic{fold}"
        or receipt["outcome_role"] != "valid"
        or receipt["mapper_outer_fold"] != fold
        or receipt["donor_test_fold"] != fold
        or receipt["genomic_test_fold"] != fold
        or receipt["valid_donor_fold"] != valid_fold
        or receipt["valid_genomic_fold"] != valid_fold
        or receipt["lineage_id"] != "hepatocyte"
        or receipt["tiles"] != 128
        or receipt["atac_channel"] != receipt["atac_reverse_complement_channel"]
        or receipt["count_like_tta_order"] != "orientation_wise_softplus_then_average"
        or (
            receipt["regional_score_names"][3]
            != "strand_tta_orientation_mean_softplus_regional_sum"
        )
        or receipt["held_ATAC_or_other_outcomes_used"]
    ):
        raise SystemExit(f"Corgi fold {fold} receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "Corgi_outcome_aligned_tile_smoke_bundle_v3",
        "model_id": "corgi_regular",
        "dataset_id": "gse296875",
        "lineage_id": "hepatocyte",
        "mapper_outer_folds": list(range(5)),
        "outcome_role": "valid",
        "crossed_split_ids": [f"donor{fold}_genomic{fold}" for fold in range(5)],
        "logical_tasks": 5,
        "outcomes_read": False,
        "open_champion_eligible": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi hepatocyte tile smoke bundle: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
