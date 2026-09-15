#!/usr/bin/env bash
#SBATCH --job-name=model-check-923
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err
set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/config-write-rollback-control-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"; FAILED="${TARGET}.failed"; ROLLBACK="${STAGING}/rollback"
cd "${ROOT}"; module purge; module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

# POSITIVE CONTROL. This job deliberately fails AFTER writing config, which is
# exactly what happened at 16:26:33 today. A rollback that never rolls back is
# a vacuous guard, so it is proven here against a scratch file rather than
# trusted because the unit tests pass.
SCRATCH="config/campaigns/rollback_control_scratch_20260825.json"
printf '{"marker": "original", "must_be_restored": true}\n' > "${ROOT}/${SCRATCH}"
ORIGINAL=$(sha256sum "${ROOT}/${SCRATCH}" | cut -d' ' -f1)
echo "scratch original digest: ${ORIGINAL}"

rollback_on_failure() {
  local rc=$?
  if [[ ${rc} -ne 0 ]]; then
    echo "run failed with ${rc}; rolling back config writes"
    if [[ -d "${ROLLBACK}" ]]; then
      PYTHONPATH=${ROOT}:${ROOT}/src python "${ROOT}/scripts/config_write_transaction.py" \
        restore --root "${ROOT}" --into "${ROLLBACK}" || echo "ROLLBACK ITSELF FAILED"
    fi
    if [[ -d "${STAGING}" && ! -e "${FAILED}" ]]; then mv "${STAGING}" "${FAILED}"; fi
  fi
}
trap rollback_on_failure EXIT

PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v tests.unit.test_config_write_transaction \
  > "${STAGING}/validation/unittest.txt" 2>&1

# snapshot BEFORE the first write
PYTHONPATH=${ROOT}:${ROOT}/src python "${ROOT}/scripts/config_write_transaction.py" \
  snapshot --root "${ROOT}" --into "${ROLLBACK}" --target "${SCRATCH}"

# write, exactly as a real writer would
printf '{"marker": "MODIFIED BY A RUN THAT WILL DIE", "must_be_restored": true}\n' > "${ROOT}/${SCRATCH}"
MODIFIED=$(sha256sum "${ROOT}/${SCRATCH}" | cut -d' ' -f1)
echo "scratch after write:     ${MODIFIED}"
[[ "${MODIFIED}" != "${ORIGINAL}" ]] || { echo "CONTROL BROKEN: the write did not change the file"; exit 3; }

echo "now failing deliberately, after the write"
exit 7
