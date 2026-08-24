#!/usr/bin/env bash
#SBATCH --job-name=masld-midas-official-probe
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
export CUBLAS_WORKSPACE_CONFIG=:4096:8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
MIDAS="${ROOT}/executions/environments/midas-v0.3.0-portable-21066061"
MIDAS_ARTIFACTS_SHA256=c861185715ed6952f8a3a09b37b04e17e5a2e8e22a23a2e3c8ab0688508fddce
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
BASE_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
TORCH_ARTIFACTS_SHA256=d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3
SITE_PACKAGES="${MIDAS}/site-packages:${TORCH}/site-packages:${BASE}/site-packages"
TARGET="${ROOT}/executions/midas-official-trainer-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed MIDAS official trainer probe: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite MIDAS official trainer probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for record in \
  "${MIDAS_ARTIFACTS_SHA256} ${MIDAS}/ARTIFACTS.json" \
  "${BASE_ARTIFACTS_SHA256} ${BASE}/ARTIFACTS.json" \
  "${TORCH_ARTIFACTS_SHA256} ${TORCH}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${MIDAS}" "${BASE}" "${TORCH}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${SITE_PACKAGES}:${ROOT}" python -m unittest -v \
  tests.unit.test_midas_official_trainer_probe \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${SITE_PACKAGES}:${ROOT}/scripts" python \
  "${ROOT}/scripts/midas_official_trainer_probe.py" \
  --output "${STAGE}/probe" --seed 8111 \
  > "${STAGE}/validation/probe.stdout.json" \
  2> "${STAGE}/validation/probe.stderr.txt"
sha256sum \
  "${ROOT}/scripts/midas_official_trainer_probe.py" \
  "${ROOT}/scripts/midas_synthetic_training_probe.py" \
  "${ROOT}/tests/unit/test_midas_official_trainer_probe.py" \
  "${ROOT}/slurm/midas_official_trainer_probe.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "probe/receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["exact_distribution"] != "scmidas 0.3.0"
    or not receipt["official_trainer_fit_passed"]
    or not receipt["safe_weights_only_state_load"]
    or not receipt["synthetic_only"]
    or receipt["project_data_read"]
    or receipt["outcomes_read"]
    or receipt["test_outcomes_read"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("MIDAS official trainer probe contract differs")
freeze_tree(root, {
    "artifact_class": "midas_official_trainer_probe",
    "model_id": "midas",
    "exact_distribution": "scmidas 0.3.0",
    "midas_runtime_artifacts_sha256": "c861185715ed6952f8a3a09b37b04e17e5a2e8e22a23a2e3c8ab0688508fddce",
    "official_trainer_fit_passed": True,
    "synthetic_only": True,
    "project_data_read": False,
    "outcomes_read": False,
    "test_outcomes_read": False,
    "champion_claim_allowed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'MIDAS official trainer probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
