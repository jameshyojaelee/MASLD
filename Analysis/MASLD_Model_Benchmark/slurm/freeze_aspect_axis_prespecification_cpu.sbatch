#!/usr/bin/env bash
#SBATCH --job-name=aspect-axis-prespec
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 0 of the four-arm aspect decomposition.
#
# Reads the four label tables, measures the tie-structure ceiling, measures the
# exact functional dependence among the labels, computes the detectable-effect
# floor, and seals the prespecification. No feature is correlated with any label
# in this stage; Stage 1 refuses to open a matrix unless this tree's digest
# matches.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
DERIVED_AFTER=${ROOT}/executions/_aspect_axis_smoke/out-21324222/directions.tsv
TARGET=${ROOT}/executions/aspect-axis-four-arm-prespec-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the prespecification.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}" "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src:${ROOT}/scripts:${ROOT}"
bash -n slurm/freeze_aspect_axis_prespecification_cpu.sbatch

"${PYTHON}" -m unittest -v tests.unit.test_aspect_axis_instrument \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" scripts/freeze_aspect_axis_prespecification.py \
  --output "${STAGE}/prespecification" \
  --derived-after-directions "${DERIVED_AFTER}" \
  > "${STAGE}/prespec_stdout.txt" 2> "${STAGE}/prespec_stderr.txt"

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum scripts/freeze_aspect_axis_prespecification.py \
  scripts/aspect_axis_common.py scripts/aspect_axis_arms.py \
  scripts/evaluate_aspect_axis_four_arm.py \
  tests/unit/test_aspect_axis_instrument.py \
  slurm/freeze_aspect_axis_prespecification_cpu.sbatch > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
import json, pathlib, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = pathlib.Path(sys.argv[1]).resolve(strict=True)
spec = json.loads((root / "prespecification" /
                   "aspect_axis_prespecification.json").read_text())
if not spec["prespec_id"].startswith("aspect_axis_count_four_arm_v"):
    raise SystemExit("wrong prespecification id")
guard = json.loads((root / "prespecification" / "ceiling_guard.json").read_text())
if not guard["passed"]:
    raise SystemExit("G_CEILING_PIN did not pass")
verify_frozen_tree(root / "prespecification")
freeze_tree(root, {"artifact_class": "aspect_axis_count_prespecification",
                   "stage": "stage_0", "status": "passed"})
verify_frozen_tree(root)
print("prespecification sealed")
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
echo "PRESPEC_DIR=${TARGET}"
sha256sum "${TARGET}/prespecification/aspect_axis_prespecification.json" \
          "${TARGET}/prespecification/ARTIFACTS.json"
