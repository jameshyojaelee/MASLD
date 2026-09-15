#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-fit-510
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${VIEWS_ROOT:?VIEWS_ROOT is required}"
: "${VIEWS_SHA256:?VIEWS_SHA256 is required}"

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
TASK=${ROOT}/config/evaluation/paired_bulk_rna_h3k27ac_task.toml
TASK_SHA256=8bf05fe693c787ffdd8fe11443e019e668ecd03318c708e390785edd837b1a28
GATE=${ROOT}/config/evaluation/paired_bulk_rna_h3k27ac_promotion_gate.json
GATE_SHA256=f26a2048346cb60f810407971c41bbc22fe862444ee5438c90da954eb44b21a1
TARGET=${ROOT}/executions/model-cpu-fit-510-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite paired baseline prediction campaign.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/predictions"
printf '%s  %s\n' "${VIEWS_SHA256}" "${VIEWS_ROOT}/ARTIFACTS.json" \
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
"${PYTHON}" - "${VIEWS_ROOT}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(root)
receipt = json.loads((root / "receipt.json").read_text())
if (
    not receipt["held_participant_h3_is_evaluator_only"]
    or receipt["outcome_or_participant_covariates_included"]
    or receipt["model_fitted"]
    or receipt["metrics_calculated"]
):
    raise SystemExit("paired baseline input view differs")
PY
for OUTER_FOLD in 0 1 2 3 4; do
  VIEW=${VIEWS_ROOT}/model/outer_${OUTER_FOLD}
  VIEW_SHA256=$(sha256sum "${VIEW}/ARTIFACTS.json" | awk '{print $1}')
  "${PYTHON}" scripts/fit_predict_gse267145_paired_baselines.py \
    --view "${VIEW}" \
    --view-artifacts-sha256 "${VIEW_SHA256}" \
    --output "${STAGE}/predictions/outer_${OUTER_FOLD}" \
    > "${STAGE}/validation/outer_${OUTER_FOLD}.stdout.json" \
    2> "${STAGE}/validation/outer_${OUTER_FOLD}.stderr.txt"
done
"${PYTHON}" - "${STAGE}/predictions" "${VIEWS_SHA256}" "${TASK_SHA256}" "${GATE_SHA256}" <<'PY'
from pathlib import Path
import csv
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive

root = Path(sys.argv[1]).resolve(strict=True)
participants = []
for fold in range(5):
    child = root / f"outer_{fold}"
    verify_frozen_tree(child)
    receipt = json.loads((child / "prediction_receipt.json").read_text())
    with (child / "query_participants.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if (
        receipt["outer_fold"] != fold
        or receipt["model_ids"] != [
            "training_mean_h3_profile",
            "pca_ridge",
            "reduced_rank_regression",
            "pls2",
            "sparse_cca",
        ]
        or receipt["held_participant_h3_read"]
        or receipt["histology_sex_age_nas_fibrosis_read"]
        or receipt["metrics_calculated"]
        or any(int(row["outer_fold"]) != fold for row in rows)
    ):
        raise SystemExit("paired fold prediction differs")
    participants.extend(row["participant_id"] for row in rows)
if len(participants) != 99 or len(set(participants)) != 99:
    raise SystemExit("paired fold prediction coverage differs")
write_json_exclusive(root / "prediction_receipt.json", {
    "schema_version": "masld-bench-gse267145-paired-baseline-campaign-v1",
    "status": "pass_development_predictions",
    "participants": 99,
    "outer_folds": 5,
    "folds_complete": True,
    "model_ids": [
        "training_mean_h3_profile",
        "pca_ridge",
        "reduced_rank_regression",
        "pls2",
        "sparse_cca",
    ],
    "input_views_artifacts_sha256": sys.argv[2],
    "task_spec_sha256": sys.argv[3],
    "promotion_gate_sha256": sys.argv[4],
    "held_participant_h3_read": False,
    "histology_sex_age_nas_fibrosis_read": False,
    "metrics_calculated": False,
    "model_selection_performed": False,
    "champion_claim_allowed": False,
})
freeze_tree(root, {
    "artifact_class": "gse267145_paired_baseline_prediction_campaign",
    "participants": 99,
    "outer_folds": 5,
    "models": 5,
    "held_participant_h3_read": False,
    "metrics_calculated": False,
    "status": "passed",
})
PY
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/build_gse267145_paired_baseline_views.py \
  scripts/fit_predict_gse267145_paired_baselines.py \
  scripts/evaluate_gse267145_paired_baseline_predictions.py \
  tests/unit/test_gse267145_paired_baseline_campaign.py \
  slurm/run_gse267145_paired_baselines_cpu.sbatch \
  "${VIEWS_ROOT}/ARTIFACTS.json" \
  "${TASK}" \
  "${GATE}" \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(root / "predictions")
freeze_tree(root, {
    "artifact_class": "gse267145_paired_baseline_fit_campaign",
    "participants": 99,
    "outer_folds": 5,
    "models": 5,
    "held_participant_h3_read": False,
    "metrics_calculated": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/predictions/ARTIFACTS.json"
