#!/usr/bin/env bash
#SBATCH --job-name=data-eval-258
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
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
FIT="${ROOT}/executions/model-fit-257-21096578"
DONOR_BIGWIGS="${ROOT}/executions/gse296875-deduplicated-tn5-bigwigs-21063992"
TRAINING_BIGWIGS="${ROOT}/executions/gse296875-training-fold-pseudobulk-21063993"
TARGET="${ROOT}/executions/data-eval-258-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi profile-head evaluation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for pair in \
  "4c46cf098c035251d0c2607e9fdb63ed703d76fb9b3b79d7982fbd32878aa353 ${FIT}/ARTIFACTS.json" \
  "01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f ${DONOR_BIGWIGS}/ARTIFACTS.json" \
  "3190aa76cd5872f50d7b0c71bcb92520f1f2428110ce2fb4c26856994bc17e5b ${TRAINING_BIGWIGS}/ARTIFACTS.json" \
  "b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 ${BASE}/ARTIFACTS.json" \
  "96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa ${IMAGE}" \
  "61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 ${CORE}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
module purge
module load singularity/4.2.2
cd "${ROOT}"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --pwd "${ROOT}" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m unittest -v tests.unit.test_evaluate_corgi_regular_crossfit_profile_head \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --pwd "${ROOT}" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python scripts/evaluate_corgi_regular_crossfit_profile_head.py \
  --fit-root "${FIT}" \
  --donor-bigwigs "${DONOR_BIGWIGS}" \
  --training-bigwigs "${TRAINING_BIGWIGS}" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/evaluation.stdout.json"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --pwd "${ROOT}" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  scripts/evaluate_corgi_regular_crossfit_profile_head.py \
  tests/unit/test_evaluate_corgi_regular_crossfit_profile_head.py \
  slurm/evaluate_corgi_regular_crossfit_profile_head_io.sbatch \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "evaluation/evaluation.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_independent_development_smoke_evaluation"
    or len(receipt["folds"]) != 5
    or receipt["histology_or_disease_labels_read"]
    or receipt["test_or_sealed_features_or_outcomes_read"]
    or receipt["promotion_gate_evaluated"]
    or receipt["full_family_ranking_performed"]
    or receipt["champion_or_external_claim_allowed"]
):
    raise SystemExit("Corgi profile-head evaluation receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "corgi_regular_crossfit_profile_head_development_evaluation",
        "model_id": "corgi_regular",
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "lineage_id": "hepatocyte",
        "adaptation_rung": "new_assay_head",
        "biological_unit": "donor",
        "genomic_block_unit": "chromosome",
        "folds": 5,
        "tests": 5,
        "development_outcomes_read": True,
        "histology_or_disease_labels_read": False,
        "test_or_sealed_features_or_outcomes_read": False,
        "promotion_gate_evaluated": False,
        "full_family_ranking_performed": False,
        "champion_or_external_claim_allowed": False,
        "development_disposition": receipt["development_disposition"],
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi cross-fit profile-head evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
