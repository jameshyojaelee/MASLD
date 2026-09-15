#!/usr/bin/env bash
#SBATCH --job-name=model-check-253
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG="${ROOT}/config/lsgkm_deltasvm_runtime_fixture.toml"
TARGET="${ROOT}/executions/model-check-253-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed CPU runtime fixture: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite CPU runtime fixture: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  "${ROOT}/scripts/acquire_lsgkm_pinned_source.py" \
  --config "${CONFIG}" \
  --output "${STAGE}/source" \
  > "${STAGE}/validation/source-acquisition.stdout.json" \
  2> "${STAGE}/validation/source-acquisition.stderr.txt"
module purge

module load GCC/13.2.0
module load zlib/1.2.13-GCCcore-13.2.0
g++ --version > "${STAGE}/validation/compiler.version.txt" 2>&1
make --version > "${STAGE}/validation/make.version.txt" 2>&1
make -C "${STAGE}/source/src" CXX=g++ gkmtrain gkmpredict \
  > "${STAGE}/validation/build.stdout.txt" \
  2> "${STAGE}/validation/build.stderr.txt"
ldd "${STAGE}/source/src/gkmtrain" \
  > "${STAGE}/validation/gkmtrain.ldd.txt" 2>&1
ldd "${STAGE}/source/src/gkmpredict" \
  > "${STAGE}/validation/gkmpredict.ldd.txt" 2>&1
sha256sum "${STAGE}/source/src/gkmtrain" "${STAGE}/source/src/gkmpredict" \
  > "${STAGE}/validation/binaries.sha256"
module purge

module load GCC/13.2.0
module load zlib/1.2.13-GCCcore-13.2.0
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
python -m unittest -v tests.unit.test_lsgkm_safe_adapter \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
python "${ROOT}/scripts/build_lsgkm_runtime_fixture.py" \
  --project-root "${ROOT}" \
  --config "${CONFIG}" \
  --source-root "${STAGE}/source" \
  --output "${STAGE}/frozen" \
  > "${STAGE}/validation/fixture.stdout.json" \
  2> "${STAGE}/validation/fixture.stderr.txt"
python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${ROOT}/config/lsgkm_deltasvm_runtime_fixture.toml" \
  "${ROOT}/scripts/acquire_lsgkm_pinned_source.py" \
  "${ROOT}/scripts/lsgkm_safe_adapter.py" \
  "${ROOT}/scripts/build_lsgkm_runtime_fixture.py" \
  "${ROOT}/tests/unit/test_lsgkm_safe_adapter.py" \
  "${ROOT}/slurm/freeze_lsgkm_runtime_fixture_cpu.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "frozen/runtime_receipt.json").read_text())
activation = json.loads(
    (root / "frozen/contract/variant_scoring_activation_plan.json").read_text()
)
if (
    receipt["status"] != "pass_cpu_runtime_and_synthetic_fixture"
    or not receipt["repeat_model_sha256_equal"]
    or not receipt["repeat_prediction_bytes_equal"]
    or not receipt["reverse_complement_fixture_pass"]
    or not receipt["direct_gkmsvm_direction_fixture_pass"]
    or not receipt["deltasvm_reverse_complement_fixture_pass"]
    or not receipt["deltasvm_sign_conversion_once"]
    or receipt["native_resume_support"] != "unsupported_retry_from_scratch_only"
    or receipt["shared_fit_readouts"] != ["direct_gkmsvm", "deltasvm"]
    or receipt["independent_model_families"] != 1
    or receipt["historical_weights_read"]
    or receipt["historical_implementation_read"]
    or receipt["outcomes_read"]
    or receipt["biological_data_read_for_training"]
    or receipt["benchmark_metrics_calculated"]
    or receipt["production_model_training_executed"]
    or receipt["production_predictions_generated"]
    or activation["status"]
    != "blocked_pending_preoutcome_roster_and_input_authority_resolution"
    or activation["authorities"]["roster_match"]
    or activation["authorities"]["crossed_outer_splits"] != 25
    or activation["provisional_shared_fit_count_if_observed_roster_is_adopted"]
    != 625
    or activation["existing_chrombpnet_negatives_admissible"]
    or activation["independent_model_families"] != 1
    or activation["production_training_executed"]
    or activation["production_predictions_generated"]
):
    raise SystemExit("LS-GKM runtime or activation contract differs")
freeze_tree(
    root,
    {
        "artifact_class": "lsgkm_cpu_runtime_and_synthetic_fixture",
        "family_id": "shared_lsgkm_gkmsvm_deltasvm",
        "source_revision": receipt["source_revision"],
        "shared_fit_readouts": ["direct_gkmsvm", "deltasvm"],
        "independent_model_families": 1,
        "synthetic_model_training_executed": True,
        "production_model_training_executed": False,
        "production_predictions_generated": False,
        "outcomes_read": False,
        "benchmark_metrics_calculated": False,
        "production_activation_status": activation["status"],
        "status": "passed_runtime_production_blocked",
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'CPU runtime fixture: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
