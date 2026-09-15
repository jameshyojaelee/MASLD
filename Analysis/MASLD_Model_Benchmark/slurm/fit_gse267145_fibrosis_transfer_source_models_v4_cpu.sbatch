#!/usr/bin/env bash
#SBATCH --job-name=model-fit-869
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
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
CONTRACT=${ROOT}/config/evaluation/gse267145_fibrosis_transfer_source_fit_v4.json
CONTRACT_SHA256=8bcf754571a8e49ba272a569d34b1b7dbae1bb9c6aade9bbf13a0b0954df2522
ACTIVATION=${ROOT}/executions/model-data-860-21129432/activation
TARGET=${ROOT}/executions/model-fit-869-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE267145 fibrosis-transfer source-fit artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CONTRACT_SHA256}" "${CONTRACT}" \
  | sha256sum --check --strict
bash -n slurm/fit_gse267145_fibrosis_transfer_source_models_v4_cpu.sbatch
# The source fit may not name, open, or hash the external expression matrix or
# the external outcome columns.  Fail loudly if a future edit adds one.
if grep -n 'gse49541_gpl570_gene_matrix\|fibrosis_stage_group\|labels\.tsv\|evaluator_only' \
    scripts/fit_gse267145_fibrosis_transfer_source_models.py "${CONTRACT}"; then
  printf 'Source fit unexpectedly names external expression or outcomes.\n' >&2
  exit 1
fi
"${PYTHON}" -m unittest -v \
  tests.unit.test_fit_gse267145_fibrosis_transfer_source_models \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/fit_gse267145_fibrosis_transfer_source_models.py \
  --benchmark-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --contract-sha256 "${CONTRACT_SHA256}" \
  --output "${STAGE}/source_fit" \
  > "${STAGE}/fit_stdout.json" 2> "${STAGE}/fit_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/fit_gse267145_fibrosis_transfer_source_models.py \
  tests/unit/test_fit_gse267145_fibrosis_transfer_source_models.py \
  slurm/fit_gse267145_fibrosis_transfer_source_models_v4_cpu.sbatch \
  "${CONTRACT}" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${ACTIVATION}/ARTIFACTS.json" \
  executions/model-data-064-21079902/fixture/molecular/ARTIFACTS.json \
  executions/model-data-061-21079623/activation/ARTIFACTS.json \
  executions/model-data-064-21079902/fixture/folds/ARTIFACTS.json \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

import numpy as np
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
fit = root / "source_fit"
receipt = json.loads((fit / "fit_receipt.json").read_text(encoding="utf-8"))
with (fit / "seed_model_index.tsv").open(encoding="utf-8", newline="") as handle:
    seeds = list(csv.DictReader(handle, delimiter="\t"))
if (
    receipt.get("status") != "passed_source_fibrosis_transfer_fit_preprocessing_locked"
    or receipt.get("fitting_participants") != 90
    or receipt.get("shared_genes") != 18147
    or receipt.get("model_seeds") != [1701, 1709, 1721, 1723, 1733]
    or receipt.get("external_labels_read") is not False
    or receipt.get("external_expression_values_read") is not False
    or receipt.get("cross_array_pooling_performed") is not False
    or receipt["endpoint_mapping_audit"]["class_counts"]
    != {"mild_f0_f1": 86, "advanced_f3_f4": 4}
    or len(seeds) != 5 * len(receipt["fitted_model_ids"])
    or sorted(receipt["fitted_model_ids"] + receipt["unfittable_model_ids"])
    != [
        "gene_median_elastic_net",
        "gene_median_linear_svm",
        "gene_median_pca_elastic_net",
        "per_array_rank_elastic_net",
    ]
    or not receipt["fitted_model_ids"]
):
    raise SystemExit("GSE267145 fibrosis source-fit verification differs")
# v1 shipped three models whose coefficients had all been shrunk to zero.  A
# constant scorer is not a baseline, so refuse to freeze one.
for row in seeds:
    if int(row["nonzero_coefficients"]) == 0 or int(row["distinct_full_fit_scores"]) < 2:
        raise SystemExit(
            f"degenerate source model reached the freeze: {row['model_id']} seed {row['model_seed']}"
        )
for row in seeds:
    state_path = fit / row["model_state_path"]
    with np.load(state_path, allow_pickle=False) as state:
        keys = set(state.files)
        if not {"selected_indices", "center", "scale"} <= keys:
            raise SystemExit("source model state lacks its preprocessing block")
        if state["selected_indices"].shape != (1000,):
            raise SystemExit("source model state selected-feature count differs")
        if row["model_id"] == "gene_median_pca_elastic_net":
            if state["pca_components"].shape != (5, 1000):
                raise SystemExit("source PCA block shape differs")
        elif "pca_components" in keys:
            raise SystemExit("a non-PCA model carries a PCA block")
        if row["model_id"] == "gene_median_linear_svm":
            if not {"svm_coef", "calibrator_coef"} <= keys:
                raise SystemExit("source SVM state lacks its calibrator")
    seed_root = state_path.parent
    freeze_tree(
        seed_root,
        {
            "artifact_class": "gse267145_fibrosis_transfer_seed_model",
            "model_id": row["model_id"],
            "model_seed": int(row["model_seed"]),
            "external_data_read": False,
            "status": "passed",
        },
    )
    verify_frozen_tree(seed_root)
freeze_tree(
    fit,
    {
        "artifact_class": "gse267145_fibrosis_transfer_source_fit",
        "fitting_participants": 90,
        "advanced_f3_f4_participants": 4,
        "shared_genes": 18147,
        "model_seeds": 5,
        "learned_models": 4,
        "selected_source_model_id": receipt["selected_source_model_id"],
        "external_expression_values_read": False,
        "external_labels_read": False,
        "status": "passed",
    },
)
verify_frozen_tree(fit)
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_fibrosis_transfer_source_fit_campaign",
        "external_expression_values_read": False,
        "external_labels_read": False,
        "external_predictions_calculated": False,
        "external_metrics_calculated": False,
        "status": "passed",
    },
)
verify_frozen_tree(root)
print(json.dumps({
    "selected_source_model_id": receipt["selected_source_model_id"],
    "training_prevalence_advanced_f3_f4": receipt["training_prevalence_advanced_f3_f4"],
    "source_model_ranking": receipt["source_model_ranking"],
    "selected_hyperparameters": receipt["selected_hyperparameters"],
}, indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/source_fit/ARTIFACTS.json"
