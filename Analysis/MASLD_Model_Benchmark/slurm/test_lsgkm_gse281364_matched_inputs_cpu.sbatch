#!/usr/bin/env bash
#SBATCH --job-name=model-check-282
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
TARGET="${ROOT}/executions/model-check-282-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite matched-input unit validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict

module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${ROOT}/src:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python -m unittest -v \
    tests.unit.test_build_lsgkm_gse281364_matched_inputs \
    > "${STAGE}/tests.stdout.txt" \
    2> "${STAGE}/tests.stderr.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:ro" \
  --env "PYTHONPATH=${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
singularity --version > "${STAGE}/singularity.version.txt"
sha256sum \
  "${ROOT}/scripts/build_lsgkm_gse281364_matched_inputs.py" \
  "${ROOT}/tests/unit/test_build_lsgkm_gse281364_matched_inputs.py" \
  "${ROOT}/slurm/test_lsgkm_gse281364_matched_inputs_cpu.sbatch" \
  > "${STAGE}/source.sha256"
module purge

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
text = (root / "tests.stderr.txt").read_text(encoding="utf-8")
if "Ran 13 tests" not in text or "OK" not in text:
    raise SystemExit("matched-input unit validation differs")
freeze_tree(root, {
    "artifact_class": "lsgkm_gse281364_matched_input_unit_validation",
    "tests": 13,
    "outcomes_read": False,
    "production_training_executed": False,
    "production_predictions_generated": False,
    "status": "passed",
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Matched-input unit validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
