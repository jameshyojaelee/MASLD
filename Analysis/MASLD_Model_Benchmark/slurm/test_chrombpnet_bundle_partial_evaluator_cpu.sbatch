#!/usr/bin/env bash
#SBATCH --job-name=model-check-286
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

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/model-check-286-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite partial evaluator tests: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

for binding in \
  "1ca00cd93df799048a6c03cb20e1136679638979fe02d4bc26a256bb2cf1453b ${ROOT}/config/evaluation/chrombpnet_bundle_21080214_partial_development.json" \
  "880b95dbd30bbe4c5697f2290efe0f4aa2e629815e72cf2e732f664b2337a42b ${ROOT}/scripts/prepare_chrombpnet_bundle_partial_evaluation.py" \
  "6770ddce7a6c068c83008b9df871c4a3fc22914a4e9df6de1bb95f90be8310f7 ${ROOT}/tests/unit/test_prepare_chrombpnet_bundle_partial_evaluation.py" \
  "702ee2bf9a02bbb5e9cb90532711bae30ba07f35a5be5d1ee11c3164419c8158 ${ROOT}/slurm/evaluate_chrombpnet_bundle_21080214_partial_cpu.sbatch"; do
  read -r expected path <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
bash -n slurm/evaluate_chrombpnet_bundle_21080214_partial_cpu.sbatch
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_prepare_chrombpnet_bundle_partial_evaluation \
  > "${STAGE}/tests.stdout.txt" \
  2> "${STAGE}/tests.stderr.txt"
python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
sha256sum \
  config/evaluation/chrombpnet_bundle_21080214_partial_development.json \
  scripts/prepare_chrombpnet_bundle_partial_evaluation.py \
  tests/unit/test_prepare_chrombpnet_bundle_partial_evaluation.py \
  slurm/evaluate_chrombpnet_bundle_21080214_partial_cpu.sbatch \
  slurm/test_chrombpnet_bundle_partial_evaluator_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
text = (root / "tests.stderr.txt").read_text(encoding="utf-8")
if "Ran 3 tests" not in text or "OK" not in text:
    raise SystemExit("partial evaluator tests differ")
freeze_tree(root, {
    "artifact_class": "chrombpnet_bundle_partial_evaluator_unit_validation",
    "tests": 3,
    "bundle_outcomes_read": False,
    "test_outcomes_read": False,
    "sealed_data_read": False,
    "status": "passed",
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Partial evaluator tests: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
