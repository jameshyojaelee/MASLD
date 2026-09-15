#!/usr/bin/env bash
#SBATCH --job-name=model-data-078
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --array=0-7%8
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.err

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
RAW_STATE_ROOT=${RAW_STATE_ROOT:-${ROOT}/executions/resumable-model-data-076/${CAMPAIGN_HASH}}
QUANT_STATE_ROOT=${ROOT}/executions/resumable-model-data-078/${CAMPAIGN_HASH}
TARGET=${ROOT}/executions/model-data-078-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
BUNDLE=$(printf '%02d' "${SLURM_ARRAY_TASK_ID}")

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite RSEM-quantification attempt artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/attempt"
install -d -m 0750 "${QUANT_STATE_ROOT}"
cd "${ROOT}"
module purge
module load RSEM/1.3.3-foss-2022a STAR/2.7.10b-GCC-11.3.0
GSE268273_TOOL_PATH=${PATH}
GSE268273_TOOL_LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}
rsem-calculate-expression --version > "${STAGE}/rsem_version.txt" 2>&1
STAR --version > "${STAGE}/star_version.txt" 2>&1
grep -Fx 'Current version: RSEM v1.3.1' "${STAGE}/rsem_version.txt"
grep -Fx '2.7.10b' "${STAGE}/star_version.txt"
module list > "${STAGE}/tool_module_list.txt" 2>&1
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
export GSE268273_TOOL_PATH
export GSE268273_TOOL_LD_LIBRARY_PATH
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
"${PYTHON}" - "${RUNTIME}" > "${STAGE}/pinned_python_runtime.json" <<'PY'
from pathlib import Path
import json
import sys
import numpy

runtime = Path(sys.argv[1]).resolve(strict=True)
numpy_path = Path(numpy.__file__).resolve(strict=True)
if not numpy_path.is_relative_to(runtime):
    raise SystemExit(f"NumPy escaped the pinned runtime: {numpy_path}")
print(json.dumps({"python": sys.version, "numpy": numpy.__version__, "numpy_path": str(numpy_path)}, sort_keys=True))
PY
bash -n slurm/quantify_gse268273_rsem_bundles_cpu.sbatch
"${PYTHON}" scripts/quantify_gse268273_rsem_bundle.py \
  --plan-root "${PLAN_ROOT}" \
  --raw-state-root "${RAW_STATE_ROOT}" \
  --reference-root "${REFERENCE_ROOT}" \
  --model-input-root "${MODEL_INPUT}" \
  --bundle-id "${SLURM_ARRAY_TASK_ID}" \
  --output "${QUANT_STATE_ROOT}/bundle_${BUNDLE}" \
  --attempt-root "${STAGE}/attempt/participants" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --seed 268273 \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
module list > "${STAGE}/python_module_list.txt" 2>&1
sha256sum \
  scripts/quantify_gse268273_rsem_bundle.py \
  slurm/quantify_gse268273_rsem_bundles_cpu.sbatch \
  "${PLAN_ROOT}/ARTIFACTS.json" \
  "${RAW_STATE_ROOT}/bundle_${BUNDLE}/ARTIFACTS.json" \
  "${REFERENCE_ROOT}/ARTIFACTS.json" \
  "${MODEL_INPUT}/ARTIFACTS.json" \
  "${QUANT_STATE_ROOT}/bundle_${BUNDLE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" "${SLURM_ARRAY_TASK_ID}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(root, {"artifact_class": "gse268273_rsem_quantification_attempt", "bundle_id": int(sys.argv[2]), "status": "passed"})
verify_frozen_tree(root)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json" "${QUANT_STATE_ROOT}/bundle_${BUNDLE}/ARTIFACTS.json"
