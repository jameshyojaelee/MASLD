#!/usr/bin/env bash
#SBATCH --job-name=model-data-888
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
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
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
ARM_ID=both_arms
PARTICIPANTS=${ROOT}/executions/model-data-072-21080651-participants
PARTICIPANTS_ARTIFACTS_SHA256=6655f23b2092e163532113e74d75d3c2a93857ee870b4f0a94fdcb4fa7f9dbc0
PREDICTIONS=${ROOT}/executions
TARGET=${ROOT}/executions/model-data-888-${SLURM_JOB_ID}-both_arms-labels
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE83452 evaluator-only label artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
bash -n slurm/build_gse83452_nash_evaluator_labels_cpu.sbatch
printf '%s  %s\n' "${PARTICIPANTS_ARTIFACTS_SHA256}" "${PARTICIPANTS}/ARTIFACTS.json" \
  | sha256sum --check --strict
install -d -m 0750 "${STAGE}/evaluator_only"
"${PYTHON}" -m unittest -v \
  tests.unit.test_build_gse83452_nash_evaluator_labels \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/build_gse83452_nash_evaluator_labels.py \
  --records "${PARTICIPANTS}/gse83452_records.tsv" \
  --prediction "${PREDICTIONS}/model-predict-886-21130315-nash_vs_nafl/predictions/training_prevalence" \
  --prediction-artifacts-sha256 108a7fde50aed871b8341c053ee0ad11b865fd649db5e77c6a1d0bdd69ffcd87 \
  --prediction "${PREDICTIONS}/model-predict-886-21130315-nash_vs_nafl/predictions/gene_median_elastic_net" \
  --prediction-artifacts-sha256 80146f0e6c1b901f453071c2497cde510940d58c2dbc15674da4649be0dab44d \
  --prediction "${PREDICTIONS}/model-predict-886-21130315-nash_vs_nafl/predictions/gene_median_linear_svm" \
  --prediction-artifacts-sha256 5bccc2e5c469f5cf8a8c7a2bbeb0849ec3ca2f84865ae7635daa51a12dc19cdd \
  --prediction "${PREDICTIONS}/model-predict-886-21130315-nash_vs_nafl/predictions/gene_median_pca_elastic_net" \
  --prediction-artifacts-sha256 0ed9bc35e8951801479f3bd802b22c76a2cc3fcca02860e8616f965152bc6b17 \
  --prediction "${PREDICTIONS}/model-predict-886-21130315-nash_vs_nafl/predictions/per_array_rank_elastic_net" \
  --prediction-artifacts-sha256 e38fdbd19810fd722ead3bb09317f29a9dcfc1676bf04c6742a7f9bd1c744156 \
  --prediction "${PREDICTIONS}/model-predict-887-21130316-nash_vs_nafl_plus_control/predictions/training_prevalence" \
  --prediction-artifacts-sha256 a2a299268ff5846b3ecba0878dcea3fef7dceff21d1ded45d97bc84998853bb6 \
  --prediction "${PREDICTIONS}/model-predict-887-21130316-nash_vs_nafl_plus_control/predictions/gene_median_elastic_net" \
  --prediction-artifacts-sha256 b0c79408cc00a3fcf9720aaabbaa664096a5d48a126143888890fbe7134514ca \
  --prediction "${PREDICTIONS}/model-predict-887-21130316-nash_vs_nafl_plus_control/predictions/gene_median_linear_svm" \
  --prediction-artifacts-sha256 9460d83c93756caf9620f1f1bc97a4b9c37e31468a4aa7b37450eb6a59b9671a \
  --prediction "${PREDICTIONS}/model-predict-887-21130316-nash_vs_nafl_plus_control/predictions/gene_median_pca_elastic_net" \
  --prediction-artifacts-sha256 ad8c82381930f36421c14efefba6a71612c5f02e72c86d9592865c8cc257dfa2 \
  --prediction "${PREDICTIONS}/model-predict-887-21130316-nash_vs_nafl_plus_control/predictions/per_array_rank_elastic_net" \
  --prediction-artifacts-sha256 802403fbbd49a8f8c9acd61c8819b3bc7766917165f563bbb863b27f5a681de6 \
  --output-labels "${STAGE}/evaluator_only/labels.tsv" \
  --output-receipt "${STAGE}/evaluator_only/label_receipt.json" \
  > "${STAGE}/label_stdout.json" 2> "${STAGE}/label_stderr.txt"
sha256sum \
  scripts/build_gse83452_nash_evaluator_labels.py \
  tests/unit/test_build_gse83452_nash_evaluator_labels.py \
  slurm/build_gse83452_nash_evaluator_labels_cpu.sbatch \
  "${PARTICIPANTS}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
evaluator = root / "evaluator_only"
receipt = json.loads((evaluator / "label_receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_evaluator_only_outcome_materialized_after_prediction_freeze"
    or receipt.get("baseline_records") != 152
    or receipt.get("endpoint_evaluable") != 148
    or receipt.get("class_counts") != {"no_nash": 44, "nash": 104}
    or receipt.get("undefined_retained_as_missing") != 4
    or receipt.get("undefined_mapped_to_no_nash") is not False
):
    raise SystemExit("GSE83452 evaluator-only label verification differs")
freeze_tree(evaluator, {
    "artifact_class": "gse83452_baseline_nash_evaluator_only_labels",
    "endpoint_evaluable": 148,
    "outcome_is_evaluator_only": True,
    "materialized_after_prediction_freeze": True,
    "status": "passed",
})
verify_frozen_tree(evaluator)
freeze_tree(root, {
    "artifact_class": "gse83452_baseline_nash_evaluator_label_campaign",
    "metrics_calculated": False,
    "status": "passed",
})
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/evaluator_only/ARTIFACTS.json"
