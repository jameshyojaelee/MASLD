#!/usr/bin/env bash
#SBATCH --job-name=data-freeze-219
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=executions/logs/data-freeze-219-%j.out
#SBATCH --error=executions/logs/data-freeze-219-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE=${ROOT}/executions/gse105127-production-9675cd48993468d4/attempts/reference-exact-21086153.failed/output.staging
TARGET=${ROOT}/executions/gse105127-reference-source-cache-${SLURM_JOB_ID}
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
cd "${ROOT}"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src:${ROOT}"
"${PYTHON}" scripts/freeze_gse105127_reference_source_cache.py \
  --source-root "${SOURCE}" --output "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
