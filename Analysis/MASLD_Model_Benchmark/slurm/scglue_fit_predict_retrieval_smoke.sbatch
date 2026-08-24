#!/usr/bin/env bash
#SBATCH --job-name=masld-scglue-retrieval-fit
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PREPARED="${ROOT}/executions/scglue-retrieval-prepared-21065511"
PREPARED_ARTIFACTS_SHA256=c79930c961cbbd44d9d198579f8b9443b6758e28cc2f4458be48d884d62a3354
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
BASE_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
GLUE="${ROOT}/executions/environments/scglue-v0.4.1-portable-21064673"
GLUE_ARTIFACTS_SHA256=371632648d40f8838f97605442623c9302886d32f84e0282f1167445ac513081
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
TORCH_ARTIFACTS_SHA256=d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3
PYTHON_SITE="${TORCH}/site-packages:${GLUE}/site-packages:${BASE}/site-packages"
TARGET="${ROOT}/executions/scglue-retrieval-fit-predict-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed scGLUE retrieval fit: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite scGLUE retrieval fit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for record in \
  "${PREPARED_ARTIFACTS_SHA256} ${PREPARED}/ARTIFACTS.json" \
  "${BASE_ARTIFACTS_SHA256} ${BASE}/ARTIFACTS.json" \
  "${GLUE_ARTIFACTS_SHA256} ${GLUE}/ARTIFACTS.json" \
  "${TORCH_ARTIFACTS_SHA256} ${TORCH}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${PREPARED}" "${BASE}" "${GLUE}" "${TORCH}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${PYTHON_SITE}:${ROOT}" python -m unittest -v \
  tests.unit.test_scglue_fit_predict_retrieval_smoke \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${PYTHON_SITE}" python \
  "${ROOT}/scripts/scglue_fit_predict_retrieval_smoke.py" \
  --inputs "${PREPARED}/prepared/model_inputs" \
  --output "${STAGE}/predictions" \
  --seed 1103 \
  > "${STAGE}/validation/fit_predict.stdout.json"
PYTHONPATH="${PYTHON_SITE}" python -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/scglue_fit_predict_retrieval_smoke.py" \
  "${ROOT}/tests/unit/test_scglue_fit_predict_retrieval_smoke.py" \
  "${ROOT}/slurm/scglue_fit_predict_retrieval_smoke.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "predictions/receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["model_roster"] != ["scglue", "paired_scglue", "linear_cca"]
    or receipt["outer_unit"] != "donor"
    or set(receipt["folds"]) != set(map(str, range(5)))
    or receipt["hidden_pair_map_read"]
    or not receipt["query_modality_ids_disjoint"]
    or receipt["retrieval_metrics_calculated"]
    or receipt["rna_to_atac_prediction_executed"]
    or receipt["sealed_rna_conditioned_atac_eligible"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("scGLUE retrieval fit receipt differs")
for fold, fold_record in receipt["folds"].items():
    if set(fold_record["models"]) != {"scglue", "paired_scglue", "linear_cca"}:
        raise SystemExit(f"scGLUE model roster differs: {fold}")
    if any(model["hidden_pair_map_read"] for model in fold_record["models"].values()):
        raise SystemExit(f"scGLUE model saw hidden pairs: {fold}")
    if any(model["retrieval_metrics_calculated"] for model in fold_record["models"].values()):
        raise SystemExit(f"scGLUE model calculated retrieval: {fold}")
freeze_tree(root, {
    "artifact_class": "scglue_retrieval_fit_predict",
    "model_roster": ["scglue", "paired_scglue", "linear_cca"],
    "dataset_id": "gse296875",
    "prepared_artifacts_sha256": "c79930c961cbbd44d9d198579f8b9443b6758e28cc2f4458be48d884d62a3354",
    "hidden_pair_map_read": False,
    "retrieval_metrics_calculated": False,
    "rna_to_atac_prediction_executed": False,
    "sealed_rna_conditioned_atac_eligible": False,
    "champion_claim_allowed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scGLUE retrieval fit/predict: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
