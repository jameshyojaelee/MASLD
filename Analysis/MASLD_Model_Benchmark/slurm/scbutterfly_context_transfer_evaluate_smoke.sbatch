#!/usr/bin/env bash
#SBATCH --job-name=masld-scbutterfly-eval
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
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
PREDICTION="${ROOT}/executions/scbutterfly-b-safe-fit-predict-21066188"
PREDICTION_ARTIFACTS_SHA256=d87dfdfcda353cf70121320d4bb428951533e4f59bd54508a175ce9e283a4033
PREPARED="${ROOT}/executions/multivi-smoke-prepared-21065093"
PREPARED_ARTIFACTS_SHA256=6c10206033545fe7ed1ecbd9f7d64e107e2b3742de06db55a838bbdb46a1cdaa
SOURCE_ROOT="${ROOT}/executions/gse296875-rna-atac-smoke-20260824T015422Z-21039865/subset"
SOURCE_ARTIFACTS_SHA256=db9e1cb847ec83a61ec79f7a206485db11e323f6b1b54fa4e10bbe4f9728aa3a
SOURCE="${SOURCE_ROOT}/gse296875_rna_atac_smoke_1000.h5"
SOURCE_SHA256=a1c05a432cec6894acae20a129ecceb153f191d9911843c86a90e13cd7f8ce14
RUNTIME="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
RUNTIME_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
SITE_PACKAGES="${RUNTIME}/site-packages"
TARGET="${ROOT}/executions/scbutterfly-b-safe-evaluation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed scButterfly evaluation: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite scButterfly evaluation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for record in \
  "${PREDICTION_ARTIFACTS_SHA256} ${PREDICTION}/ARTIFACTS.json" \
  "${PREPARED_ARTIFACTS_SHA256} ${PREPARED}/ARTIFACTS.json" \
  "${SOURCE_ARTIFACTS_SHA256} ${SOURCE_ROOT}/ARTIFACTS.json" \
  "${SOURCE_SHA256} ${SOURCE}" \
  "${RUNTIME_ARTIFACTS_SHA256} ${RUNTIME}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${PREDICTION}" "${PREPARED}" "${SOURCE_ROOT}" "${RUNTIME}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${SITE_PACKAGES}:${ROOT}" python -m unittest -v \
  tests.unit.test_evaluate_multivi_smoke \
  tests.unit.test_evaluate_scbutterfly_context_transfer_smoke \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${SITE_PACKAGES}" python \
  "${ROOT}/scripts/evaluate_scbutterfly_context_transfer_smoke.py" \
  --prediction-root "${PREDICTION}/predictions" \
  --prepared "${PREPARED}/prepared" \
  --source-h5 "${SOURCE}" \
  --source-h5-sha256 "${SOURCE_SHA256}" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/evaluation.stdout.json"
sha256sum \
  "${ROOT}/scripts/evaluate_scbutterfly_context_transfer_smoke.py" \
  "${ROOT}/scripts/evaluate_multivi_smoke.py" \
  "${ROOT}/tests/unit/test_evaluate_multivi_smoke.py" \
  "${ROOT}/tests/unit/test_evaluate_scbutterfly_context_transfer_smoke.py" \
  "${ROOT}/slurm/scbutterfly_context_transfer_evaluate_smoke.sbatch" \
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
    or receipt["model_id"] != "scbutterfly_b_safe"
    or receipt["evaluation_role"] != "development_smoke"
    or receipt["pairing_topology"] != "same_nucleus_training"
    or receipt["biological_unit"] != "donor"
    or receipt["n_donors"] != 39
    or receipt["n_lineages"] != 5
    or receipt["n_regions"] != 10000
    or not receipt["prediction_frozen_before_outcomes"]
    or receipt["held_atac_exposed_to_model"]
    or not receipt["held_atac_read_by_evaluator_only"]
    or receipt["test_outcomes_read"]
    or receipt["external_evaluation"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("scButterfly evaluation contract differs")
freeze_tree(root, {
    "artifact_class": "scbutterfly_b_safe_development_evaluation",
    "model_id": "scbutterfly_b_safe",
    "dataset_id": "gse296875",
    "prediction_artifacts_sha256": "d87dfdfcda353cf70121320d4bb428951533e4f59bd54508a175ce9e283a4033",
    "prediction_frozen_before_outcomes": True,
    "held_atac_exposed_to_model": False,
    "test_outcomes_read": False,
    "rna_conditioned_atac_eligible": True,
    "sealed_inference_eligible": False,
    "champion_claim_allowed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scButterfly-B safe development evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
