#!/usr/bin/env bash
#SBATCH --job-name=model-activation-081
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
CONTRACT=${ROOT}/config/evaluation/gse260666_external_histology_activation_v1.json
CONTRACT_SHA256=26d66a94df0a2fd9c6f8632f6500bc5db74e4ea23ae919bd0d824c73aeea96f5
TARGET=${ROOT}/executions/model-activation-081-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE260666 activation artifacts.\n' >&2
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
bash -n slurm/freeze_gse260666_external_histology_activation_cpu.sbatch
"${PYTHON}" -m unittest -v \
  tests.unit.test_audit_gse260666_external_histology_activation \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/audit_gse260666_external_histology_activation.py \
  --benchmark-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --contract-sha256 "${CONTRACT_SHA256}" \
  --output "${STAGE}/activation" \
  > "${STAGE}/audit_stdout.json" 2> "${STAGE}/audit_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/audit_gse260666_external_histology_activation.py \
  tests/unit/test_audit_gse260666_external_histology_activation.py \
  slurm/freeze_gse260666_external_histology_activation_cpu.sbatch \
  "${CONTRACT}" \
  "${RUNTIME}/ARTIFACTS.json" \
  executions/model-select-080-21092216/selection/ARTIFACTS.json \
  executions/model-training-069-production-v4/ARTIFACTS.json \
  executions/model-data-064-21079902/fixture/molecular/ARTIFACTS.json \
  executions/model-data-068-21080339/ARTIFACTS.json \
  executions/model-data-068-21080339/fixture/model_input/ARTIFACTS.json \
  config/evaluation/gse260666_external_histology_transfer_task.toml \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
activation = root / "activation"
receipt = json.loads((activation / "activation_receipt.json").read_text(encoding="utf-8"))
with (activation / "common_stable_gene_axis.tsv").open(encoding="utf-8", newline="") as handle:
    common = list(csv.DictReader(handle, delimiter="\t"))
with (activation / "source_stage3_recipe_distribution.tsv").open(encoding="utf-8", newline="") as handle:
    recipes = list(csv.DictReader(handle, delimiter="\t"))
if (
    receipt.get("status") != "passed_donor_safe_external_development_activation"
    or receipt.get("selected_source_model_id") != "rna_hvg_pca_linear_svm"
    or receipt.get("common_stable_genes") != 20410
    or receipt.get("external_labels_read") is not False
    or receipt.get("endpoint_mapping_registered_before_source_selection") is not True
    or receipt.get("project_sealed") is not False
    or len(common) != 20410
    or sum(row["locked_full_fit_recipe"] == "true" for row in recipes) != 1
):
    raise SystemExit("GSE260666 activation verification differs")
freeze_tree(
    activation,
    {
        "artifact_class": "gse260666_external_histology_activation",
        "selected_source_model_id": "rna_hvg_pca_linear_svm",
        "common_stable_genes": 20410,
        "external_labels_read": False,
        "project_sealed": False,
        "status": "passed",
    },
)
verify_frozen_tree(activation)
freeze_tree(
    root,
    {
        "artifact_class": "gse260666_external_histology_activation_campaign",
        "external_feature_values_read": False,
        "external_labels_read": False,
        "model_fitted": False,
        "external_predictions_calculated": False,
        "status": "passed",
    },
)
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/activation/ARTIFACTS.json"
