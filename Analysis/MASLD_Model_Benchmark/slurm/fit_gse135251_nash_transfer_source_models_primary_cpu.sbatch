#!/usr/bin/env bash
#SBATCH --job-name=model-fit-882
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
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
CONTRACT=${ROOT}/config/evaluation/gse135251_nash_transfer_source_fit_nash_vs_nafl_v1.json
CONTRACT_SHA256=a6c6cf3bcdf1dbc123344f85be227e54c5a3094a216fb05a71bb1704e65c9687
SOURCE=${ROOT}/executions/model-data-880-21130257-gse135251-source/source
ACTIVATION=${ROOT}/executions/model-data-881-21130281-activation/activation
TARGET=${ROOT}/executions/model-fit-882-${SLURM_JOB_ID}-${ARM_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE135251 source-fit artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CONTRACT_SHA256}" "${CONTRACT}" | sha256sum --check --strict
bash -n slurm/fit_gse135251_nash_transfer_source_models_primary_cpu.sbatch
# The source fit must never name the external expression matrix, the external
# record table, or the evaluator outcome.
if grep -n 'gse83452_gpl16686_gene_matrix\|gse83452_records\|baseline_participant_masks\|labels\.tsv' \
    scripts/fit_gse135251_nash_transfer_source_models.py "${CONTRACT}"; then
  printf 'The source fit unexpectedly names a GSE83452 outcome-bearing artifact.\n' >&2
  exit 1
fi
"${PYTHON}" -m unittest -v \
  tests.unit.test_fit_gse135251_nash_transfer_source_models \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/fit_gse135251_nash_transfer_source_models.py \
  --benchmark-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --contract-sha256 "${CONTRACT_SHA256}" \
  --output "${STAGE}/source_fit" \
  > "${STAGE}/fit_stdout.json" 2> "${STAGE}/fit_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/fit_gse135251_nash_transfer_source_models.py \
  tests/unit/test_fit_gse135251_nash_transfer_source_models.py \
  slurm/fit_gse135251_nash_transfer_source_models_primary_cpu.sbatch \
  "${CONTRACT}" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE}/ARTIFACTS.json" \
  "${ACTIVATION}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" "${ARM_ID}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
arm_id = sys.argv[2]
fit = root / "source_fit"
receipt = json.loads((fit / "fit_receipt.json").read_text(encoding="utf-8"))
expected = {"no_nash": 41, "nash": 131} if arm_id == "nash_vs_nafl" else {"no_nash": 49, "nash": 131}
if (
    receipt.get("status") != "passed_source_nash_transfer_fit_preprocessing_locked"
    or receipt.get("arm_id") != arm_id
    or receipt.get("shared_genes") != 21465
    or receipt.get("endpoint_mapping_audit", {}).get("class_counts") != expected
    or receipt.get("external_labels_read") is not False
    or receipt.get("external_expression_values_read") is not False
    or receipt.get("cross_array_pooling_performed") is not False
    or receipt.get("age_sex_logistic_fitted") is not False
    or "age_sex_logistic" not in receipt.get("unfittable_model_ids", [])
    or receipt.get("selected_source_model_id") in receipt.get("unfittable_model_ids", [])
):
    raise SystemExit("GSE135251 source-fit verification differs")
for model_id in receipt["fitted_model_ids"]:
    for seed in receipt["model_seeds"]:
        seed_root = fit / "models" / model_id / f"seed_{seed}"
        seed_receipt = json.loads((seed_root / "receipt.json").read_text(encoding="utf-8"))
        if seed_receipt["nonzero_coefficients"] <= 0 or seed_receipt["distinct_full_fit_scores"] < 2:
            raise SystemExit(f"{model_id} seed {seed} is degenerate")
        freeze_tree(
            seed_root,
            {
                "artifact_class": "gse135251_nash_transfer_source_model_state",
                "arm_id": arm_id,
                "model_id": model_id,
                "model_seed": seed,
                "external_labels_read": False,
                "status": "passed",
            },
        )
        verify_frozen_tree(seed_root)
freeze_tree(
    fit,
    {
        "artifact_class": "gse135251_nash_transfer_source_fit",
        "arm_id": arm_id,
        "fitting_participants": receipt["fitting_participants"],
        "shared_genes": receipt["shared_genes"],
        "selected_source_model_id": receipt["selected_source_model_id"],
        "external_labels_read": False,
        "external_expression_values_read": False,
        "source_roster_admission_state": "proposed_not_approved",
        "status": "passed",
    },
)
verify_frozen_tree(fit)
freeze_tree(
    root,
    {
        "artifact_class": "gse135251_nash_transfer_source_fit_campaign",
        "arm_id": arm_id,
        "external_labels_read": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
verify_frozen_tree(root)
print(json.dumps({
    "arm_id": receipt["arm_id"],
    "fitting_participants": receipt["fitting_participants"],
    "training_prevalence_nash": receipt["training_prevalence_nash"],
    "fitted_model_ids": receipt["fitted_model_ids"],
    "unfittable_model_ids": receipt["unfittable_model_ids"],
    "selected_source_model_id": receipt["selected_source_model_id"],
    "source_model_ranking": receipt["source_model_ranking"],
    "no_source_model_exceeds_its_own_null": receipt["no_source_model_exceeds_its_own_null"],
    "training_partition_linearly_separable": receipt["training_partition_linearly_separable"],
}, indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/source_fit/ARTIFACTS.json"
