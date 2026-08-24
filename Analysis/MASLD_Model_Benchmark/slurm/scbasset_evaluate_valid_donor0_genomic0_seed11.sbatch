#!/usr/bin/env bash
#SBATCH --job-name=masld-scbasset-d0g0-s11-eval
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PREDICTION="${ROOT}/executions/scbasset-donor0-genomic0-seed11-valid-21064859"
PREDICTION_SHA=c27f080243968ec5815e158fc22280c1718ef9b0cb0ea08012442755d35a9db0
INPUT="${ROOT}/executions/scbasset-donor0-genomic0-inputs-21064455"
INPUT_SHA=02f76f4f7af62063a606943e2b1bcf11ca0c8d8a63f29879ed31bfc1fbc30dc0
BIGWIGS="${ROOT}/executions/gse296875-deduplicated-tn5-bigwigs-21063992"
BIGWIGS_SHA=01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
SPLIT_SHA=00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
TARGET="${ROOT}/executions/scbasset-donor0-genomic0-seed11-valid-evaluation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" ]]; then
    if [[ -e "${failed}" ]]; then
      printf 'Cannot preserve failed scBasset evaluation: %s\n' "${failed}" >&2
      return
    fi
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed scBasset evaluation: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite scBasset evaluation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for record in \
  "${PREDICTION_SHA} ${PREDICTION}/ARTIFACTS.json" \
  "${INPUT_SHA} ${INPUT}/ARTIFACTS.json" \
  "${BIGWIGS_SHA} ${BIGWIGS}/ARTIFACTS.json" \
  "${SPLIT_SHA} ${SPLIT}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${PREDICTION}" "${INPUT}" "${BIGWIGS}" "${SPLIT}" <<'PY'
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
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python -m unittest -v tests.unit.test_evaluate_scbasset_valid \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python "${ROOT}/scripts/evaluate_scbasset_valid.py" \
    --prediction-root "${PREDICTION}" \
    --prediction-artifacts-sha256 "${PREDICTION_SHA}" \
    --input-root "${INPUT}" \
    --input-artifacts-sha256 "${INPUT_SHA}" \
    --donor-bigwigs "${BIGWIGS}" \
    --donor-bigwigs-artifacts-sha256 "${BIGWIGS_SHA}" \
    --split-root "${SPLIT}" \
    --split-artifacts-sha256 "${SPLIT_SHA}" \
    --output "${STAGE}/evaluation" \
    --workers 8 \
    > "${STAGE}/validation/evaluation.stdout.json"
sha256sum \
  "${ROOT}/scripts/evaluate_scbasset_valid.py" \
  "${ROOT}/tests/unit/test_evaluate_scbasset_valid.py" \
  "${ROOT}/slurm/scbasset_evaluate_valid_donor0_genomic0_seed11.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "evaluation/evaluation.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["evaluation_role"] != "valid"
    or receipt["n_donors"] != 9
    or receipt["n_lineages"] != 5
    or receipt["n_regions"] != 16000
    or receipt["held_valid_atac_exposed_to_model"]
    or receipt["test_regions_read"]
    or receipt["test_atac_read"]
    or receipt["test_predictions_read"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("scBasset validation evaluation contract differs")
freeze_tree(root, {
    "artifact_class": "scbasset_valid_development_evaluation",
    "model_id": "scbasset",
    "dataset_id": "gse296875",
    "split_id": "donor0_genomic0",
    "seed": 11,
    "evaluation_role": "valid",
    "prediction_artifacts_sha256": receipt["prediction_artifacts_sha256"],
    "test_outcomes_read": False,
    "champion_claim_allowed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scBasset validation evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
