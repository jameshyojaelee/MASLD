#!/usr/bin/env bash
#SBATCH --job-name=corgi-regional-eval
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${PREDICTIONS_ROOT:?PREDICTIONS_ROOT is required}"
: "${PREDICTIONS_SHA256:?PREDICTIONS_SHA256 is required}"
[[ "${PREDICTIONS_ROOT}" =~ /executions/corgi-hepatocyte-tile-smoke-bundle-v3-[0-9]+$ ]] || {
  printf 'Corgi prediction root contract differs\n' >&2; exit 1;
}
[[ "${PREDICTIONS_SHA256}" =~ ^[0-9a-f]{64}$ ]] || {
  printf 'Corgi prediction manifest hash differs\n' >&2; exit 1;
}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
DONOR_BIGWIGS="${ROOT}/executions/gse296875-deduplicated-tn5-bigwigs-21063992"
TRAINING_BIGWIGS="${ROOT}/executions/gse296875-training-fold-pseudobulk-21063993"
UNIT_TEST="${ROOT}/executions/corgi-regional-evaluator-test-21070309"
TARGET="${ROOT}/executions/corgi-regional-development-smoke-evaluation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi evaluation artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for pair in \
  "${PREDICTIONS_SHA256} ${PREDICTIONS_ROOT}/ARTIFACTS.json" \
  "b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 ${BASE}/ARTIFACTS.json" \
  "96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa ${IMAGE}" \
  "61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 ${CORE}/ARTIFACTS.json" \
  "01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f ${DONOR_BIGWIGS}/ARTIFACTS.json" \
  "3190aa76cd5872f50d7b0c71bcb92520f1f2428110ce2fb4c26856994bc17e5b ${TRAINING_BIGWIGS}/ARTIFACTS.json" \
  "47a17200eb098796f2fde6976e02927fc3fb68103cd1eb7906453eb144bbb26b ${UNIT_TEST}/ARTIFACTS.json" \
  "8008d5a8b8803b61890db4823d4b27be4c472f44881ca3788ce728f01be8f731 ${ROOT}/scripts/evaluate_corgi_regional_smoke.py" \
  "f2725c66fae08cb4e2002008441a36e30b8591abc6cdc11731fff1fa3893fa16 ${ROOT}/tests/unit/test_evaluate_corgi_regional_smoke.py"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${PREDICTIONS_ROOT}" "${UNIT_TEST}" "${BASE}" "${CORE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]: verify_frozen_tree(Path(value).resolve(strict=True))
PY
module purge
module load singularity/4.2.2
cd "${ROOT}"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --pwd "${ROOT}" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m unittest -v tests.unit.test_evaluate_corgi_regional_smoke \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --pwd "${ROOT}" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python scripts/evaluate_corgi_regional_smoke.py \
  --predictions "${PREDICTIONS_ROOT}" \
  --donor-bigwigs "${DONOR_BIGWIGS}" \
  --training-bigwigs "${TRAINING_BIGWIGS}" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/evaluation.stdout.json"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --pwd "${ROOT}" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
value = json.loads((root / "evaluation/evaluation.json").read_text())
if (
    value["schema_version"] != "masld-bench-corgi-regional-smoke-evaluation-v1"
    or value["status"] != "pass_development_smoke_evaluation"
    or len(value["folds"]) != 5
    or len(value["context_arm_summaries"]) != 5
    or not value["development_outcomes_read"]
    or not value["held_validation_atac_read"]
    or value["held_test_or_sealed_outcomes_read"]
    or value["histology_or_disease_labels_read"]
    or value["champion_claim_allowed"]
    or value["promotion_gate_evaluated"]
):
    raise SystemExit("Corgi development smoke evaluation contract differs")
PY
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --pwd "${ROOT}" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  scripts/evaluate_corgi_regional_smoke.py \
  tests/unit/test_evaluate_corgi_regional_smoke.py \
  slurm/evaluate_corgi_regional_smoke.sbatch \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${PREDICTIONS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
value = json.loads((root / "evaluation/evaluation.json").read_text())
freeze_tree(root, {
    "artifact_class": "corgi_regional_development_smoke_evaluation",
    "dataset_id": "gse296875", "model_id": "corgi_regular",
    "lineage_id": "hepatocyte", "prediction_artifacts_sha256": sys.argv[2],
    "biological_unit": "donor", "genomic_block_unit": "chromosome",
    "development_outcomes_read": True, "held_test_or_sealed_outcomes_read": False,
    "histology_or_disease_labels_read": False, "champion_claim_allowed": False,
    "promotion_gate_evaluated": False, "logical_units": sum(
        item["units"] for item in value["context_arm_summaries"].values()
    ), "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi regional development smoke evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
