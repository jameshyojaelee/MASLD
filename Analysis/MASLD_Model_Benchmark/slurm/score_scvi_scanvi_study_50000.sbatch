#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-score-503
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
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${PREDICTIONS_ROOT:?PREDICTIONS_ROOT must be the frozen ensemble directory}"
: "${PREDICTIONS_SHA256:?PREDICTIONS_SHA256 must bind its ARTIFACTS.json}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
TARGET=${ROOT}/executions/model-cpu-score-503-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite scVI/scANVI evaluation artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${SOURCE_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SPLIT_SHA256}" "${SPLIT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${PREDICTIONS_SHA256}" "${PREDICTIONS_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict

cd "${ROOT}"
export PYTHONPATH=${ROOT}:${ROOT}/src
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${PYTHON}" -m unittest \
  tests.unit.test_score_cell_baselines_study_50000 \
  tests.unit.test_aggregate_scvi_scanvi_study_50000 \
  > "${STAGE}/validation/unit-tests.txt" \
  2>&1
"${PYTHON}" scripts/score_scvi_scanvi_study_50000.py \
  --source "${SOURCE}" \
  --split "${SPLIT}" \
  --predictions "${PREDICTIONS_ROOT}" \
  --output "${STAGE}/evaluation" \
  --expected-source-artifacts-sha256 "${SOURCE_SHA256}" \
  --expected-split-artifacts-sha256 "${SPLIT_SHA256}" \
  --expected-prediction-artifacts-sha256 "${PREDICTIONS_SHA256}" \
  > "${STAGE}/validation/evaluator.stdout.json" \
  2> "${STAGE}/validation/evaluator.stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  scripts/score_scvi_scanvi_study_50000.py \
  scripts/score_cell_baselines_study_50000.py \
  tests/unit/test_score_cell_baselines_study_50000.py \
  tests/unit/test_aggregate_scvi_scanvi_study_50000.py \
  slurm/score_scvi_scanvi_study_50000.sbatch \
  "${SOURCE}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  "${PREDICTIONS_ROOT}/ARTIFACTS.json" \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

# Two decisions a reader of these numbers has to have. Written here in the
# wrapper rather than inside the scorer, because the scorer is bound by SHA-256
# and adding a field to it would move that pin. freeze_tree below covers this
# file, so it is part of the sealed record rather than a note beside it.
cat > "${STAGE}/recorded_decisions.json" <<'DECISIONS'
{
  "schema_version": "masld-bench-scvi-scanvi-evaluation-recorded-decisions-v1",
  "why_this_file_exists": "These are decisions taken with their costs stated, not rules that went missing. A reader who finds them only in a message thread cannot audit them.",

  "hardware_comparability_exception": {
    "decision": "scVI and scANVI were fit on B6K. They are compared with the five HVG/PCA baselines, which were fit on CPU. This is a recorded exception to the letter of config/b6k_routing_guidance_20260825.json, taken deliberately on 2026-08-30.",
    "the_rule": "That record says a B6K set may only be compared with a B6K set, on the same pool and one frozen runtime.",
    "why_the_exception_holds": "The rule exists to stop numerical drift across CUDA runtime and GPU architecture. The five HVG/PCA baselines are sklearn under the frozen cpu_baseline_smoke environment lock, with no CUDA in the path at all, so they carry no exposure to the mechanism the rule guards.",
    "what_enforcing_the_letter_would_have_cost": "Re-running those five on the CPUs of a B6K node. Roughly 8.5 hours; it would hold one of only 22 B6K GPUs idle to run sklearn, because B6K nodes sit in the gpu partition alone; and it would very likely NOT reproduce the published numbers, because the originals ran on ne1dc4-023 which is Intel v4 while B6K nodes are AMD. Threaded BLAS changes reduction order across that boundary, and this cluster has already been bitten by it once. A non-reproducing re-run would replace rather than confirm the baselines, moving the strongest baseline for a reason that is not any model improving.",
    "conclusion": "The compliance was more expensive and less valid than the exception.",
    "what_this_exception_does_NOT_cover": "TF-Sapiens. It was fit on L40S under torch 2.5.1 with cu124, verified from jobs 21080428 and 21082720 and the runtime pip freeze. cu124 cannot execute on B6K at all, so that arm could never be reproduced there. This is a feasibility exclusion and no exception is available for it."
  },

  "celltypist_baseline_is_not_designed": {
    "state": "NOT DESIGNED. Not merely unrun.",
    "evidence": "config/models/mandatory_baselines.toml declares celltypist_baseline with implementation_type train_locally. No scored output exists anywhere under executions.",
    "still_to_be_frozen_before_it_could_run": ["normalization", "feature selection", "solver", "regularization", "class balancing", "seed", "majority-vote policy", "calibration"],
    "consequence": "cell_state_mapping names eight baselines. Five are fit on CPU, two are these scVI and scANVI results on B6K, and the eighth has no design. The set is structurally incomplete for a reason no amount of compute fixes, and independent of any hardware question."
  },

  "reading_these_numbers_against_tf_sapiens": {
    "rule": "If scVI or scANVI lands above 0.91887, the strongest baseline moves and TF-Sapiens two_layer_mlp at 0.93887 no longer clears its required 0.02 margin over it.",
    "incumbent_strongest_baseline": "hvg_pca_knn at 0.9074773235260757, from executions/model-cpu-score-502-21097344/evaluation/strongest_baseline.json",
    "requirement": "Recomputing a candidate's margin against a changed baseline set is a separate deliberate act. When it happens, report both margins side by side with the baseline named in each. A verdict that changes because a baseline improved is legitimate; one that changes silently is not.",
    "seed_axis": "These scVI and scANVI results and the five HVG/PCA baselines use seeds 20260824, 20260825 and 20260826. TF-Sapiens used 1103, 1201 and 1301. That is a further axis on which the two sets are not interchangeable, independent of both silicon and runtime."
  }
}
DECISIONS

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(root / "evaluation")
receipt = json.loads(
    (root / "evaluation/evaluation_receipt.json").read_text(encoding="utf-8")
)
if (
    receipt.get("status") != "pass_development_evaluation_only"
    or receipt.get("models") != ["scvi_baseline", "scanvi_baseline"]
    or receipt.get("rows") != 50_000
    or receipt.get("donors") != 102
    or receipt.get("studies") != 7
    or receipt.get("models_fit") is not False
    or receipt.get("preprocessing_fit") is not False
    or receipt.get("calibration_fit") is not False
    or receipt.get("bootstrap", {}).get("resamples") != 10_000
    or receipt.get("bootstrap", {}).get("cells_are_independent_units") is not False
    or receipt.get("histology_read") is not False
    or receipt.get("sealed_outcomes_read") is not False
):
    raise SystemExit("scVI/scANVI evaluation receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "scvi_scanvi_study_heldout_evaluation_campaign",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "split_id": "resource_atlas_study_outer_5fold_v1",
        "models_fit": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/evaluation/ARTIFACTS.json"
