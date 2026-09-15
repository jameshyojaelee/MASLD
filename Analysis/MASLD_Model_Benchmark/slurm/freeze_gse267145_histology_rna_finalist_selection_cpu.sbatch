#!/usr/bin/env bash
#SBATCH --job-name=model-select-080
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
CONTRACT=${ROOT}/config/evaluation/gse267145_histology_rna_finalist_selection_v1.json
CONTRACT_SHA256=bb80cac22cfa0c4142a3f89d4865c7196a21d8522c7d4549bda831e7634c4fcd
TARGET=${ROOT}/executions/model-select-080-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite RNA finalist selection artifacts.\n' >&2
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
bash -n slurm/freeze_gse267145_histology_rna_finalist_selection_cpu.sbatch
"${PYTHON}" -m unittest -v \
  tests.unit.test_audit_gse267145_histology_rna_finalist_selection \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/audit_gse267145_histology_rna_finalist_selection.py \
  --benchmark-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --contract-sha256 "${CONTRACT_SHA256}" \
  --output "${STAGE}/selection" \
  > "${STAGE}/audit_stdout.json" 2> "${STAGE}/audit_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/audit_gse267145_histology_rna_finalist_selection.py \
  tests/unit/test_audit_gse267145_histology_rna_finalist_selection.py \
  slurm/freeze_gse267145_histology_rna_finalist_selection_cpu.sbatch \
  "${CONTRACT}" \
  "${RUNTIME}/ARTIFACTS.json" \
  executions/model-training-069-production-v4/ARTIFACTS.json \
  executions/model-scoring-078-21092130/ARTIFACTS.json \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
selection = root / "selection"
receipt = json.loads((selection / "selection_receipt.json").read_text(encoding="utf-8"))
with (selection / "rna_candidate_primary_calibration.tsv").open(encoding="utf-8", newline="") as handle:
    candidates = list(csv.DictReader(handle, delimiter="\t"))
with (selection / "stage3_hyperparameter_one_se_audit.tsv").open(encoding="utf-8", newline="") as handle:
    one_se = list(csv.DictReader(handle, delimiter="\t"))
if (
    receipt.get("status") != "passed_frozen_development_rna_finalist_selection"
    or receipt.get("selected_model_id") != "rna_hvg_pca_linear_svm"
    or receipt.get("secondary_endpoint_values_used_for_selection") is not False
    or receipt.get("champion_claim_allowed") is not False
    or len(candidates) != 4
    or sum(row["selected_rna_finalist"] == "true" for row in candidates) != 1
    or len(one_se) != 100
    or any(row["rederived_match"] != "true" for row in one_se)
):
    raise SystemExit("RNA finalist selection verification differs")
freeze_tree(
    selection,
    {
        "artifact_class": "gse267145_histology_rna_finalist_selection",
        "selected_model_id": "rna_hvg_pca_linear_svm",
        "candidate_models": 4,
        "stage3_one_standard_error_audits": 100,
        "secondary_endpoints_used_for_selection": False,
        "development_only": True,
        "status": "passed",
    },
)
verify_frozen_tree(selection)
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_histology_rna_finalist_selection_campaign",
        "selected_model_id": "rna_hvg_pca_linear_svm",
        "external_data_accessed": False,
        "model_fitted": False,
        "external_predictions_calculated": False,
        "status": "passed",
    },
)
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/selection/ARTIFACTS.json"
