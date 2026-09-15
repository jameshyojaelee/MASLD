#!/usr/bin/env bash
#SBATCH --job-name=model-eval-626
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
SCRIPT="${ROOT}/scripts/evaluate_lsgkm_gse281364_dinucleotide_scores.py"
OUTCOME_HELPER="${ROOT}/scripts/evaluate_gse281364_dna_lm_common_lane.py"
TEST="${ROOT}/tests/unit/test_evaluate_lsgkm_gse281364_dinucleotide_scores.py"
LOCK="${ROOT}/executions/model-check-625-21100532/ARTIFACTS.json"
LOCK_SHA256=de04d91ea4b790e3de8f6b8bab2509d2b0ca3e9311f901f32982250f8e512d1f
TARGET="${ROOT}/executions/model-eval-626-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite LS-GKM development evaluation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

printf '%s  %s\n' f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb "${RUNTIME}" | sha256sum --check --strict
printf '%s  %s\n' e9678f9c4c2919a65d67167f441b3938c5ad972b29b870269989969cfd8c0343 "${CONFIG}" | sha256sum --check --strict
printf '%s  %s\n' 69d366cf002ab59e041e2b1eb9f33bbe9e3059ac56c91bdbf4983c35647e886d "${SCRIPT}" | sha256sum --check --strict
printf '%s  %s\n' 39730c63b7370ea05153482f881da5165287e0a8ff1ad757a1d9399db2f7ab82 "${OUTCOME_HELPER}" | sha256sum --check --strict
printf '%s  %s\n' 437377bd633418582edfc307328e89cd8af8f2b0c672bd0c036322f63b667883 "${TEST}" | sha256sum --check --strict
printf '%s  %s\n' "${LOCK_SHA256}" "${LOCK}" | sha256sum --check --strict

cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_evaluate_lsgkm_gse281364_dinucleotide_scores \
  > "${STAGE}/tests.stdout.txt" 2> "${STAGE}/tests.stderr.txt"
"${PYTHON}" "${SCRIPT}" \
  --project-root "${ROOT}" --config "${CONFIG}" \
  --selection-lock-artifacts "${LOCK}" \
  --selection-lock-artifacts-sha256 "${LOCK_SHA256}" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/evaluation.stdout.json" 2> "${STAGE}/evaluation.stderr.txt"
"${PYTHON}" -m pip freeze > "${STAGE}/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/python.version.txt" 2>&1
sha256sum "${SCRIPT}" "${OUTCOME_HELPER}" "${TEST}" "${CONFIG}" "${LOCK}" "${RUNTIME}" \
  "${ROOT}/slurm/evaluate_lsgkm_gse281364_dinucleotide_scores_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "evaluation/evaluation_receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_locked_exposed_development_evaluation"
    or receipt["candidate_count"] != 3
    or receipt["prediction_rows"] != 30990
    or receipt["long_range_blocks"] != 239
    or len(receipt["results"]) != 3
    or receipt["readouts_are_independent_model_families"]
    or receipt["mandatory_baselines_complete"]
    or receipt["candidate_selection_or_shortlist_created"]
    or receipt["external_evaluation"]
    or receipt["champion_claim"]
    or receipt["confirmatory_p_values_reported"]
    or receipt["sealed_assets_read"]
    or receipt["outcomes_read_in_isolated_evaluator"] is not True
    or receipt["control_predictions_read_in_isolated_evaluator"] is not True
):
    raise SystemExit("LS-GKM development evaluation receipt differs")
freeze_tree(root, {
    "artifact_class": "lsgkm_gse281364_exact_dinucleotide_development_evaluation",
    "dataset_id": "gse281364",
    "task_id": "variant_to_regulation",
    "candidate_count": receipt["candidate_count"],
    "prediction_rows": receipt["prediction_rows"],
    "elements": receipt["elements"],
    "long_range_blocks": receipt["long_range_blocks"],
    "fixed_seeds": receipt["fixed_seeds"],
    "readouts_are_independent_model_families": False,
    "mandatory_baselines_complete": False,
    "candidate_selection_or_shortlist_created": False,
    "external_evaluation": False,
    "champion_claim": False,
    "sealed_assets_read": False,
    "status": receipt["status"],
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'LS-GKM development evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
