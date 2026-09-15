#!/usr/bin/env bash
#SBATCH --job-name=model-fit-082
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
CONTRACT=${ROOT}/config/evaluation/gse267145_rna_full_fit_for_gse260666_v1.json
CONTRACT_SHA256=9c2753bc3e64fb0b0c96f4623a8a302a26865392cb586c4c073b4750f6465730
TARGET=${ROOT}/executions/model-fit-082-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE267145 RNA source full-fit artifacts.\n' >&2
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
bash -n slurm/fit_gse267145_rna_full_fit_for_gse260666_cpu.sbatch
"${PYTHON}" -m unittest -v \
  tests.unit.test_fit_gse267145_rna_full_fit_for_gse260666 \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/fit_gse267145_rna_full_fit_for_gse260666.py \
  --benchmark-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --contract-sha256 "${CONTRACT_SHA256}" \
  --output "${STAGE}/source_fit" \
  > "${STAGE}/fit_stdout.json" 2> "${STAGE}/fit_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/fit_gse267145_rna_full_fit_for_gse260666.py \
  tests/unit/test_fit_gse267145_rna_full_fit_for_gse260666.py \
  slurm/fit_gse267145_rna_full_fit_for_gse260666_cpu.sbatch \
  "${CONTRACT}" \
  "${RUNTIME}/ARTIFACTS.json" \
  executions/model-activation-081-21092266/activation/ARTIFACTS.json \
  executions/model-select-080-21092216/selection/ARTIFACTS.json \
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
with (fit / "selected_feature_axis.tsv").open(encoding="utf-8", newline="") as handle:
    features = list(csv.DictReader(handle, delimiter="\t"))
with (fit / "seed_model_index.tsv").open(encoding="utf-8", newline="") as handle:
    seeds = list(csv.DictReader(handle, delimiter="\t"))
for row in seeds:
    state = np.load(fit / row["model_state_path"], allow_pickle=False)
    if (
        state["center"].shape != (1000,)
        or state["pca_components"].shape != (5, 1000)
        or state["svm_coef"].shape != (3, 5)
        or state["svm_calibrator_coef"].shape != (3, 3)
        or state["elastic_coef"].shape != (3, 5)
        or state["centroid_centroids"].shape != (3, 5)
    ):
        raise SystemExit("source full-fit model-state shape differs")
if (
    receipt.get("status") != "passed_source_full99_fit_preprocessing_locked"
    or receipt.get("participants") != 99
    or receipt.get("model_seeds") != [1701, 1709, 1721, 1723, 1733]
    or receipt.get("external_expression_values_read") is not False
    or receipt.get("external_labels_read") is not False
    or len(features) != 1000
    or len(seeds) != 5
):
    raise SystemExit("source full-fit verification differs")
for row in seeds:
    seed_root = fit / Path(row["model_state_path"]).parent
    freeze_tree(
        seed_root,
        {
            "artifact_class": "gse267145_rna_full_fit_seed_model",
            "model_seed": int(row["model_seed"]),
            "external_data_read": False,
            "status": "passed",
        },
    )
    verify_frozen_tree(seed_root)
freeze_tree(
    fit,
    {
        "artifact_class": "gse267145_rna_full_fit_for_gse260666",
        "participants": 99,
        "model_seeds": 5,
        "selected_model_id": "rna_hvg_pca_linear_svm",
        "context_models": 3,
        "external_expression_values_read": False,
        "external_labels_read": False,
        "status": "passed",
    },
)
verify_frozen_tree(fit)
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_rna_full_fit_campaign",
        "external_expression_values_read": False,
        "external_labels_read": False,
        "external_predictions_calculated": False,
        "external_metrics_calculated": False,
        "status": "passed",
    },
)
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/source_fit/ARTIFACTS.json"
