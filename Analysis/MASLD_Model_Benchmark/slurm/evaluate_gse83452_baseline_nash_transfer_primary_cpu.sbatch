#!/usr/bin/env bash
#SBATCH --job-name=model-eval-890
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
ARM_ID=nash_vs_nafl
GATE=${ROOT}/config/evaluation/gse83452_baseline_nash_transfer_promotion_gate.json
GATE_SHA256=f8f28bf4812fd124e198b0bf810cbafd2775ccb9e86bb1f301a717a76e17e1e7
EVALUATOR_ONLY=${ROOT}/executions/model-data-888-21130317-both_arms-labels/evaluator_only
EVALUATOR_ONLY_ARTIFACTS_SHA256=a325746d4188794392bd25e545f664bad3be5581409c13b16bb050bca78287e8
PREDICTIONS=${ROOT}/executions/model-predict-886-21130315-nash_vs_nafl/predictions
TARGET=${ROOT}/executions/model-eval-890-${SLURM_JOB_ID}-nash_vs_nafl-eval
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE83452 evaluation artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
bash -n slurm/evaluate_gse83452_baseline_nash_transfer_primary_cpu.sbatch
printf '%s  %s\n' "${GATE_SHA256}" "${GATE}" | sha256sum --check --strict
"${PYTHON}" -m unittest -v \
  tests.unit.test_evaluate_gse83452_baseline_nash_transfer \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/evaluate_gse83452_baseline_nash_transfer.py \
  --evaluator-only "${EVALUATOR_ONLY}" \
  --evaluator-only-artifacts-sha256 "${EVALUATOR_ONLY_ARTIFACTS_SHA256}" \
  --gate "${GATE}" \
  --gate-sha256 "${GATE_SHA256}" \
  --candidate-model-id gene_median_elastic_net \
  --arm-id nash_vs_nafl \
  --gate-eligible-arm true \
  --unfittable-baseline age_sex_logistic \
  --prediction "${PREDICTIONS}/training_prevalence" \
  --prediction-artifacts-sha256 108a7fde50aed871b8341c053ee0ad11b865fd649db5e77c6a1d0bdd69ffcd87 \
  --prediction "${PREDICTIONS}/gene_median_elastic_net" \
  --prediction-artifacts-sha256 80146f0e6c1b901f453071c2497cde510940d58c2dbc15674da4649be0dab44d \
  --prediction "${PREDICTIONS}/gene_median_linear_svm" \
  --prediction-artifacts-sha256 5bccc2e5c469f5cf8a8c7a2bbeb0849ec3ca2f84865ae7635daa51a12dc19cdd \
  --prediction "${PREDICTIONS}/gene_median_pca_elastic_net" \
  --prediction-artifacts-sha256 0ed9bc35e8951801479f3bd802b22c76a2cc3fcca02860e8616f965152bc6b17 \
  --prediction "${PREDICTIONS}/per_array_rank_elastic_net" \
  --prediction-artifacts-sha256 e38fdbd19810fd722ead3bb09317f29a9dcfc1676bf04c6742a7f9bd1c744156 \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/evaluation_stdout.json" 2> "${STAGE}/evaluation_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/evaluate_gse83452_baseline_nash_transfer.py \
  tests/unit/test_evaluate_gse83452_baseline_nash_transfer.py \
  slurm/evaluate_gse83452_baseline_nash_transfer_primary_cpu.sbatch \
  "${GATE}" \
  "${EVALUATOR_ONLY}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" "${ARM_ID}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
arm_id = sys.argv[2]
evaluation = root / "evaluation"
receipt = json.loads((evaluation / "receipt.json").read_text(encoding="utf-8"))
verdict = json.loads((evaluation / "promotion_gate_verdict.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "passed_external_development_scoring"
    or receipt.get("arm_id") != arm_id
    or receipt.get("participants") != 148
    or receipt.get("frozen_predictions") != 152
    or receipt.get("undefined_mapped_to_no_nash") is not False
    or receipt.get("prevalence_is_not_the_auprc_null") is not True
    or receipt.get("champion_claim_eligible") is not False
    or verdict.get("verdict") not in {"PASS", "FAIL"}
):
    raise SystemExit("GSE83452 evaluation verification differs")
freeze_tree(evaluation, {
    "artifact_class": "gse83452_baseline_nash_transfer_evaluation",
    "arm_id": arm_id,
    "participants": 148,
    "promotion_gate_verdict": verdict["verdict"],
    "external_development_only": True,
    "champion_eligible": False,
    "status": "passed",
})
verify_frozen_tree(evaluation)
freeze_tree(root, {
    "artifact_class": "gse83452_baseline_nash_transfer_evaluation_campaign",
    "arm_id": arm_id,
    "promotion_gate_verdict": verdict["verdict"],
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps(verdict, indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/evaluation/ARTIFACTS.json"
cat "${TARGET}/evaluation/promotion_gate_verdict.json"
