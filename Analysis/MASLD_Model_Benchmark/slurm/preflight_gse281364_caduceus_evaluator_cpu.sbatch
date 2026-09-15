#!/usr/bin/env bash
#SBATCH --job-name=evaluator-preflight
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
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
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON="${PROJECT}/.agent/scanpy_env/bin/python"
TARGET="${ROOT}/executions/caduceus-evaluator-preflight-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Caduceus evaluator preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for binding in \
  "bb837a6abbf6d01a504d9539345c970d1fec14d37f19b6a7b75431251dd3e034 ${ROOT}/scripts/evaluate_gse281364_caduceus_common_lane.py" \
  "77d9bdfeedde8ad55ca5858fc7325cdb80fe59e1d27d9737abe8eb32b34d1b7c ${ROOT}/tests/unit/test_evaluate_gse281364_caduceus_common_lane.py" \
  "39730c63b7370ea05153482f881da5165287e0a8ff1ad757a1d9399db2f7ab82 ${ROOT}/scripts/evaluate_gse281364_dna_lm_common_lane.py" \
  "066d37b0497008caac7b96006ec2b2ebab132afa06189bb518cab435a1ab9afc ${ROOT}/tests/unit/test_evaluate_gse281364_dna_lm_common_lane.py"; do
  read -r expected path <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
cd "${ROOT}"
PYTHONPATH="${ROOT}" "${PYTHON}" -m unittest -v \
  tests.unit.test_evaluate_gse281364_dna_lm_common_lane \
  tests.unit.test_evaluate_gse281364_caduceus_common_lane \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}" "${PYTHON}" scripts/evaluate_gse281364_caduceus_common_lane.py --help \
  > "${STAGE}/validation/cli-help.stdout.txt" \
  2> "${STAGE}/validation/cli-help.stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  scripts/evaluate_gse281364_caduceus_common_lane.py \
  tests/unit/test_evaluate_gse281364_caduceus_common_lane.py \
  scripts/evaluate_gse281364_dna_lm_common_lane.py \
  tests/unit/test_evaluate_gse281364_dna_lm_common_lane.py \
  slurm/preflight_gse281364_caduceus_evaluator_cpu.sbatch \
  > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
tests = (root / "validation/tests.stderr.txt").read_text()
help_text = (root / "validation/cli-help.stdout.txt").read_text()
if (
    "Ran 9 tests" not in tests
    or "OK" not in tests
    or "--embeddings" not in help_text
    or "--outcomes" not in help_text
    or (root / "validation/cli-help.stderr.txt").stat().st_size != 0
):
    raise SystemExit("Caduceus evaluator preflight differs")
freeze_tree(root, {
    "artifact_class": "gse281364_caduceus_evaluator_outcome_free_preflight",
    "dataset_id": "gse281364",
    "model_id": "caduceus",
    "synthetic_tests": 9,
    "project_outcomes_read": False,
    "prediction_artifact_read": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Caduceus evaluator preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
