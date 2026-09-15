#!/usr/bin/env bash
#SBATCH --job-name=data-check-218
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=executions/logs/data-check-218-%j.out
#SBATCH --error=executions/logs/data-check-218-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CAMPAIGN=${ROOT}/executions/gse105127-production-9675cd48993468d4
TARGET=${ROOT}/executions/data-check-218-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
cd "${ROOT}"
if [[ -e ${TARGET} || -e ${STAGE} ]]; then
  printf 'Refusing to overwrite RNA remote identity census.\n' >&2
  exit 1
fi
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src:${ROOT}"
"${PYTHON}" -m unittest -v tests.unit.test_diagnose_gse105127_rna_remote_identity
"${PYTHON}" scripts/diagnose_gse105127_rna_remote_identity.py \
  --plan-root "${CAMPAIGN}/plan" \
  --local-root "${CAMPAIGN}/rna_raw" \
  --output "${STAGE}"
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
