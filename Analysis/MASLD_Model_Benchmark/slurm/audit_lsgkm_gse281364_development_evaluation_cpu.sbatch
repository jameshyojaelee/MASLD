#!/usr/bin/env bash
#SBATCH --job-name=model-check-627
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
RUNTIME="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/ARTIFACTS.json"
CONFIG="${ROOT}/config/lsgkm_gse281364_dinucleotide_development_evaluation.json"
SCRIPT="${ROOT}/scripts/audit_lsgkm_gse281364_development_evaluation.py"
TEST="${ROOT}/tests/unit/test_audit_lsgkm_gse281364_development_evaluation.py"
EVALUATION="${ROOT}/executions/model-eval-626-21100570/ARTIFACTS.json"
EVALUATION_SHA256=3a2db575d98c1fb22954c9f993ff3ba6390204437f5d285d0fda5ff5aa54276e
TARGET="${ROOT}/executions/model-check-627-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite LS-GKM evaluation audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

printf '%s  %s\n' f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb "${RUNTIME}" | sha256sum --check --strict
printf '%s  %s\n' e9678f9c4c2919a65d67167f441b3938c5ad972b29b870269989969cfd8c0343 "${CONFIG}" | sha256sum --check --strict
printf '%s  %s\n' 9940ffc7ec20c1ba39622d79b1c9180d48ebadfc292b0727c53c2d0adf0e430e "${SCRIPT}" | sha256sum --check --strict
printf '%s  %s\n' d33c6f93b1422b6b90800b1e0c5ba1454d40b4e09e9fea90476a234227a92117 "${TEST}" | sha256sum --check --strict
printf '%s  %s\n' "${EVALUATION_SHA256}" "${EVALUATION}" | sha256sum --check --strict

cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v tests.unit.test_audit_lsgkm_gse281364_development_evaluation \
  > "${STAGE}/tests.stdout.txt" 2> "${STAGE}/tests.stderr.txt"
"${PYTHON}" "${SCRIPT}" \
  --project-root "${ROOT}" --config "${CONFIG}" \
  --evaluation-artifacts "${EVALUATION}" \
  --evaluation-artifacts-sha256 "${EVALUATION_SHA256}" \
  --output "${STAGE}/audit" \
  > "${STAGE}/audit.stdout.json" 2> "${STAGE}/audit.stderr.txt"
"${PYTHON}" -m pip freeze > "${STAGE}/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/python.version.txt" 2>&1
sha256sum "${SCRIPT}" "${TEST}" "${CONFIG}" "${EVALUATION}" "${RUNTIME}" \
  "${ROOT}/slurm/audit_lsgkm_gse281364_development_evaluation_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "audit/evaluation_audit.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_independent_exposed_development_evaluation_rederivation"
    or receipt["candidate_count"] != 3
    or receipt["prediction_rows"] != 30990
    or receipt["elements"] != 1033
    or receipt["long_range_blocks"] != 239
    or receipt["bootstrap_resamples"] != 10000
    or len(receipt["results"]) != 3
    or receipt["candidate_selection_or_shortlist_created"]
    or receipt["external_evaluation"]
    or receipt["champion_claim"]
    or receipt["confirmatory_p_values_reported"]
    or receipt["sealed_assets_read"]
):
    raise SystemExit("LS-GKM evaluation audit receipt differs")
freeze_tree(root, {
    "artifact_class": "lsgkm_gse281364_exposed_development_evaluation_audit",
    "evaluation_artifacts_sha256": receipt["evaluation_artifacts_sha256"],
    "candidate_count": receipt["candidate_count"],
    "prediction_rows": receipt["prediction_rows"],
    "elements": receipt["elements"],
    "long_range_blocks": receipt["long_range_blocks"],
    "bootstrap_resamples": receipt["bootstrap_resamples"],
    "candidate_selection_or_shortlist_created": False,
    "external_evaluation": False,
    "champion_claim": False,
    "sealed_assets_read": False,
    "status": receipt["status"],
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'LS-GKM evaluation audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
