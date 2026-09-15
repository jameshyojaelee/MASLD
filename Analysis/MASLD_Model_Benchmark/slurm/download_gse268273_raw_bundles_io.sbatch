#!/usr/bin/env bash
#SBATCH --job-name=model-data-076
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --array=0-7%8
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${PLAN_ROOT:?PLAN_ROOT must name the frozen model-data-075 plan directory}"

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
CAMPAIGN_HASH=$("${PYTHON}" -c 'import json,sys; print(json.load(open(sys.argv[1]))["campaign_hash"])' "${PLAN_ROOT}/receipt.json")
RAW_STATE_ROOT=${ROOT}/executions/resumable-model-data-076/${CAMPAIGN_HASH}
TARGET=${ROOT}/executions/model-data-076-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite raw-download attempt artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
install -d -m 0750 "${RAW_STATE_ROOT}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
bash -n slurm/download_gse268273_raw_bundles_io.sbatch
"${PYTHON}" scripts/download_gse268273_raw_bundle.py \
  --plan-root "${PLAN_ROOT}" \
  --state-root "${RAW_STATE_ROOT}" \
  --bundle-id "${SLURM_ARRAY_TASK_ID}" \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/download_gse268273_raw_bundle.py \
  slurm/download_gse268273_raw_bundles_io.sbatch \
  "${PLAN_ROOT}/ARTIFACTS.json" \
  "${RAW_STATE_ROOT}/bundle_$(printf '%02d' "${SLURM_ARRAY_TASK_ID}")/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" "${SLURM_ARRAY_TASK_ID}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(root, {"artifact_class": "gse268273_raw_download_attempt", "bundle_id": int(sys.argv[2]), "status": "passed"})
verify_frozen_tree(root)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json"
