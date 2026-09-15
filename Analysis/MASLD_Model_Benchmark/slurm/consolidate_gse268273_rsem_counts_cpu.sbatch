#!/usr/bin/env bash
#SBATCH --job-name=model-data-079
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${PLAN_ROOT:?PLAN_ROOT must name the frozen model-data-075 plan directory}"
: "${REFERENCE_ROOT:?REFERENCE_ROOT must name the frozen model-data-077 artifact}"

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
MODEL_INPUT=${ROOT}/executions/model-data-074-21080832/fixture/model_input
CAMPAIGN_HASH=$("${PYTHON}" -c 'import json,sys; print(json.load(open(sys.argv[1]))["campaign_hash"])' "${PLAN_ROOT}/receipt.json")
QUANT_STATE_ROOT=${QUANT_STATE_ROOT:-${ROOT}/executions/resumable-model-data-078/${CAMPAIGN_HASH}}
TARGET=${ROOT}/executions/model-data-079-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite consolidated RSEM-count artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
bash -n slurm/consolidate_gse268273_rsem_counts_cpu.sbatch
"${PYTHON}" scripts/consolidate_gse268273_rsem_counts.py \
  --plan-root "${PLAN_ROOT}" \
  --quant-state-root "${QUANT_STATE_ROOT}" \
  --model-input-root "${MODEL_INPUT}" \
  --reference-root "${REFERENCE_ROOT}" \
  --bundles 8 \
  --output "${STAGE}/model_input" \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/consolidate_gse268273_rsem_counts.py \
  slurm/consolidate_gse268273_rsem_counts_cpu.sbatch \
  "${PLAN_ROOT}/ARTIFACTS.json" \
  "${REFERENCE_ROOT}/ARTIFACTS.json" \
  "${MODEL_INPUT}/ARTIFACTS.json" \
  "${STAGE}/model_input/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "model_input/receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "passed_outcome_free_quantification_training_transform_blocked"
    or receipt.get("participants") != 109
    or receipt.get("target_stable_genes") != 14078
    or receipt.get("target_duplicate_sum_groups") != 11
    or receipt.get("labels_included") is not False
    or receipt.get("fit_or_score_performed") is not False
):
    raise SystemExit("consolidated GSE268273 RSEM receipt differs")
freeze_tree(root, {"artifact_class": "gse268273_raw_rsem_consolidation_execution", "status": "passed"})
verify_frozen_tree(root)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/model_input/ARTIFACTS.json"
