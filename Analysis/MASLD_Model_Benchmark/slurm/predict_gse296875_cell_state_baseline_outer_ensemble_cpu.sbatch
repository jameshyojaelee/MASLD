#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-transfer-336
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
MICROMAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
RUNTIME_PREFIX=/gpfs/commons/home/jameslee/micromamba/envs/immune-atlas
FEATURES=${ROOT}/executions/model-data-332-21100842/view/features
FEATURES_SHA256=1cacb52032beb27969c69ae2d00cfdf3d18b333f1458bb4fb830f5225d79dc6c
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
BASELINE=${ROOT}/executions/model-cpu-train-503-21097161/predictions
TARGET=${ROOT}/executions/model-cpu-transfer-336-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
trap 'status=$?; if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi; exit ${status}' EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then exit 1; fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
module purge
module load python/3.11.5-GCCcore-13.2.0
python - <<'PY' "${FEATURES}" "${SOURCE}" "${SPLIT}" "${BASELINE}"
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]: verify_frozen_tree(Path(value))
PY
BASELINE_SHA256=$(sha256sum "${BASELINE}/ARTIFACTS.json" | cut -d ' ' -f 1)
for pair in "${FEATURES_SHA256} ${FEATURES}/ARTIFACTS.json" "${SOURCE_SHA256} ${SOURCE}/ARTIFACTS.json" "${SPLIT_SHA256} ${SPLIT}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
module purge
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env \
  PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}:${ROOT}/src" \
  OMP_NUM_THREADS=${OMP_NUM_THREADS} MKL_NUM_THREADS=${MKL_NUM_THREADS} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS} \
  python scripts/predict_gse296875_cell_state_baseline_outer_ensemble.py \
    --features "${FEATURES}" --expected-features-sha256 "${FEATURES_SHA256}" \
    --source "${SOURCE}" --expected-source-sha256 "${SOURCE_SHA256}" \
    --split "${SPLIT}" --expected-split-sha256 "${SPLIT_SHA256}" \
    --baseline "${BASELINE}" --expected-baseline-sha256 "${BASELINE_SHA256}" \
    --output "${STAGE}/branch" \
  > "${STAGE}/validation/predict.stdout.json" 2> "${STAGE}/validation/predict.stderr.txt"
"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env PYTHONNOUSERSITE=1 \
  python -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum scripts/predict_gse296875_cell_state_baseline_outer_ensemble.py \
  slurm/predict_gse296875_cell_state_baseline_outer_ensemble_cpu.sbatch \
  "${FEATURES}/ARTIFACTS.json" "${SOURCE}/ARTIFACTS.json" "${SPLIT}/ARTIFACTS.json" \
  "${BASELINE}/ARTIFACTS.json" > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
python - <<'PY' "${STAGE}"
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]); meta = verify_frozen_tree(root / "branch")["metadata"]
if meta.get("branch_id") != "task_native_outer_ensemble" or meta.get("models") is None or len(meta["models"]) != 5:
    raise SystemExit("baseline transfer branch differs")
freeze_tree(root, {
    "artifact_class": "gse296875_cell_state_task_native_transfer_campaign",
    "view_id": "gse296875_rna_cell_state_external_development_7500_v1",
    "models": meta["models"],
    "rows": 7500,
    "query_labels_read": False,
    "query_donor_ids_read": False,
    "query_atac_read": False,
    "metrics_calculated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/branch/ARTIFACTS.json"
