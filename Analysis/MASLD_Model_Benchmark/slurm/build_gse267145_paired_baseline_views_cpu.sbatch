#!/usr/bin/env bash
#SBATCH --job-name=model-data-087
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
MOLECULAR=${ROOT}/executions/model-data-064-21079902/fixture/molecular
MOLECULAR_SHA256=2a60ef5d5d904f5f833833ebef75e7946e7db3b6b9aaf502cb84ddedf97eb456
FOLDS=${ROOT}/executions/model-data-064-21079902/fixture/folds
FOLDS_SHA256=9f5b96e69f217ba643b4b2fd4f3e165049b0cc2bfe9d65c99a7789816b4d7ee3
TASK=${ROOT}/config/evaluation/paired_bulk_rna_h3k27ac_task.toml
TASK_SHA256=8bf05fe693c787ffdd8fe11443e019e668ecd03318c708e390785edd837b1a28
GATE=${ROOT}/config/evaluation/paired_bulk_rna_h3k27ac_promotion_gate.json
GATE_SHA256=f26a2048346cb60f810407971c41bbc22fe862444ee5438c90da954eb44b21a1
TARGET=${ROOT}/executions/model-data-087-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite paired baseline views.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${MOLECULAR_SHA256}" "${MOLECULAR}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${FOLDS_SHA256}" "${FOLDS}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${TASK_SHA256}" "${TASK}" | sha256sum --check --strict
printf '%s  %s\n' "${GATE_SHA256}" "${GATE}" | sha256sum --check --strict
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict

cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${PYTHON}" -m unittest -v tests.unit.test_gse267145_paired_baseline_campaign \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" scripts/build_gse267145_paired_baseline_views.py \
  --molecular "${MOLECULAR}" \
  --molecular-artifacts-sha256 "${MOLECULAR_SHA256}" \
  --folds "${FOLDS}" \
  --folds-artifacts-sha256 "${FOLDS_SHA256}" \
  --task-spec "${TASK}" \
  --task-spec-sha256 "${TASK_SHA256}" \
  --promotion-gate "${GATE}" \
  --promotion-gate-sha256 "${GATE_SHA256}" \
  --output "${STAGE}/views" \
  > "${STAGE}/validation/builder_stdout.json" \
  2> "${STAGE}/validation/builder_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/build_gse267145_paired_baseline_views.py \
  scripts/fit_predict_gse267145_paired_baselines.py \
  scripts/evaluate_gse267145_paired_baseline_predictions.py \
  tests/unit/test_gse267145_paired_baseline_campaign.py \
  slurm/build_gse267145_paired_baseline_views_cpu.sbatch \
  slurm/run_gse267145_paired_baselines_cpu.sbatch \
  slurm/evaluate_gse267145_paired_baselines_cpu.sbatch \
  "${MOLECULAR}/ARTIFACTS.json" \
  "${FOLDS}/ARTIFACTS.json" \
  "${TASK}" \
  "${GATE}" \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(root / "views")
receipt = json.loads((root / "views/receipt.json").read_text())
if (
    receipt["participants"] != 99
    or receipt["outer_folds"] != 5
    or receipt["pairing"] != "same_sample_different_aliquot"
    or receipt["rna_features"] != 42_163
    or receipt["h3k27ac_opaque_features"] != 96_460
    or receipt["masked_retired_rna_genes"] != 1_122
    or not receipt["held_participant_h3_is_evaluator_only"]
    or receipt["outcome_or_participant_covariates_included"]
    or receipt["coordinate_semantics_inferred"]
    or receipt["model_fitted"]
    or receipt["metrics_calculated"]
):
    raise SystemExit("paired baseline view receipt differs")
freeze_tree(root, {
    "artifact_class": "gse267145_paired_baseline_view_campaign",
    "participants": 99,
    "outer_folds": 5,
    "held_participant_h3_is_evaluator_only": True,
    "model_fitted": False,
    "metrics_calculated": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum \
  "${TARGET}/ARTIFACTS.json" \
  "${TARGET}/views/ARTIFACTS.json" \
  "${TARGET}/views/model/ARTIFACTS.json" \
  "${TARGET}/views/evaluator/ARTIFACTS.json"
