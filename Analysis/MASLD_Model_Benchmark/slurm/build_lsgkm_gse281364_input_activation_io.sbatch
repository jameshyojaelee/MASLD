#!/usr/bin/env bash
#SBATCH --job-name=model-data-283
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=48:00:00
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
UNIT="${ROOT}/executions/model-check-282-21097187"
CONFIG="${ROOT}/config/lsgkm_gse281364_input_activation.json"
REFERENCE_ROOT=/gpfs/commons/home/jameslee/reference_genome
TARGET="${ROOT}/executions/model-data-283-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite LS-GKM input activation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 4536034aa0d99e2b4b58b520fa92b04d36e645c600e89c43c0774521b0fc15e9 \
  "${UNIT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 5b9e93dddd6a0f32b112fa8096945d65bf601753c21613c70663385dfc2228f3 \
  "${ROOT}/scripts/build_lsgkm_gse281364_matched_inputs.py" | sha256sum --check --strict
printf '%s  %s\n' 44d0f9d9f20392d768f8d37b1881298c17660aea667e9e9ca34e3fc7834e2fab \
  "${ROOT}/scripts/audit_lsgkm_gse281364_input_activation.py" | sha256sum --check --strict
printf '%s  %s\n' 89b7240acdfc66c47f5ac61be0139215d9b584b1203650e009765e571b02411a \
  "${ROOT}/tests/unit/test_build_lsgkm_gse281364_matched_inputs.py" | sha256sum --check --strict
printf '%s  %s\n' eac6951cd4763bd0a13ec3e562d843c61025f174ea6b0703e501ef8d37aa071c \
  "${CONFIG}" | sha256sum --check --strict

module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --bind "${REFERENCE_ROOT}:${REFERENCE_ROOT}:ro" \
  --env "PYTHONPATH=${ROOT}:${ROOT}/src:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python -c \
    'from pathlib import Path; import sys; compile(Path(sys.argv[1]).read_text(encoding="utf-8"), sys.argv[1], "exec")' \
    "${ROOT}/scripts/build_lsgkm_gse281364_matched_inputs.py"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --bind "${REFERENCE_ROOT}:${REFERENCE_ROOT}:ro" \
  --env "PYTHONPATH=${ROOT}:${ROOT}/src:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python "${ROOT}/scripts/build_lsgkm_gse281364_matched_inputs.py" \
    --project-root "${ROOT}" \
    --config "${CONFIG}" \
    --output "${STAGE}/activation" \
    > "${STAGE}/builder.stdout.json" \
    2> "${STAGE}/builder.stderr.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:ro" \
  --env "PYTHONPATH=${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
singularity --version > "${STAGE}/singularity.version.txt"
module purge

module load python/3.11.5-GCCcore-13.2.0
python -c \
  'from pathlib import Path; import sys; compile(Path(sys.argv[1]).read_text(encoding="utf-8"), sys.argv[1], "exec")' \
  "${ROOT}/scripts/audit_lsgkm_gse281364_input_activation.py"
python "${ROOT}/scripts/audit_lsgkm_gse281364_input_activation.py" \
  --activation "${STAGE}/activation" \
  --output "${STAGE}/independent_audit.json" \
  > "${STAGE}/audit.stdout.json" \
  2> "${STAGE}/audit.stderr.txt"
sha256sum \
  "${ROOT}/scripts/build_lsgkm_gse281364_matched_inputs.py" \
  "${ROOT}/scripts/audit_lsgkm_gse281364_input_activation.py" \
  "${ROOT}/tests/unit/test_build_lsgkm_gse281364_matched_inputs.py" \
  "${ROOT}/slurm/build_lsgkm_gse281364_input_activation_io.sbatch" \
  "${CONFIG}" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
contract = json.loads(
    (root / "activation/contract/production_activation_contract.json").read_text(
        encoding="utf-8"
    )
)
audit = json.loads((root / "independent_audit.json").read_text(encoding="utf-8"))
if (
    contract["status"] != "inputs_frozen_production_fit_not_executed"
    or contract["split_design"]["shared_fits"] != 25
    or audit["status"] != "pass_outcome_blind_input_activation"
    or audit["fit_count"] != 25
    or audit["matched_pairs_minimum"] < 10_000
    or contract["outcomes_read"]
    or contract["reporter_counts_read"]
    or contract["sealed_assets_read"]
    or contract["production_training_executed"]
    or contract["production_predictions_generated"]
):
    raise SystemExit("LS-GKM input activation audit differs")
freeze_tree(root, {
    "artifact_class": "lsgkm_gse281364_outcome_blind_production_input_activation",
    "shared_fit_count": 25,
    "active_state_roster": ["hepatocyte"],
    "outcomes_read": False,
    "reporter_counts_read": False,
    "sealed_assets_read": False,
    "production_training_executed": False,
    "production_predictions_generated": False,
    "status": "passed_inputs_frozen_production_fit_not_executed",
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'LS-GKM input activation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
