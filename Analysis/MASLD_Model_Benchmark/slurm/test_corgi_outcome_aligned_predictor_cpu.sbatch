#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-predict-test
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/corgi-nvidia-pytorch-24.01-py3-recovery-21064690"
IMAGE="${RUNTIME}/nvidia-pytorch-24.01-py3.sif"
OVERLAY="${ROOT}/executions/environments/corgi-rotary-emb-0.1-flashattn2.0.4-build-r1-21064845"
TARGET="${ROOT}/executions/corgi-outcome-aligned-predictor-test-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi predictor test artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in \
  "01de778c17c67d92de9a2e242350e17443dbfe76b7c5a5f8bc447eae5408cd3c ${RUNTIME}/ARTIFACTS.json" \
  "c60c4ebf0b32df7d158a3a1d781f5d45fdf4187ab0fba9f4bea35605a92e4c3a ${OVERLAY}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
install -d -m 0750 "${STAGE}/validation"
module load singularity/4.2.2
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  --env "PYTHONPATH=${ROOT}:${OVERLAY}/site-packages" \
  "${IMAGE}" python -m unittest \
    tests.unit.test_corgi_predict_outcome_aligned_tiles \
    > "${STAGE}/validation/tests.stdout.txt" \
    2> "${STAGE}/validation/tests.stderr.txt"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/corgi_predict_outcome_aligned_tiles.py" \
  "${ROOT}/scripts/corgi_predict_valid_contexts.py" \
  "${ROOT}/scripts/corgi_probe_native_runtime.py" \
  "${ROOT}/tests/unit/test_corgi_predict_outcome_aligned_tiles.py" \
  "${ROOT}/slurm/test_corgi_outcome_aligned_predictor_cpu.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
test_log = (root / "validation/tests.stderr.txt").read_text()
if "\nOK\n" not in test_log or "FAILED" in test_log or "ERROR:" in test_log:
    raise SystemExit("Corgi outcome-aligned predictor unit test did not pass")
freeze_tree(root, {
    "artifact_class": "corgi_outcome_aligned_predictor_cpu_test",
    "model_id": "corgi_regular",
    "outcomes_read": False,
    "gpu_used": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi outcome-aligned predictor CPU test: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
