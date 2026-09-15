#!/usr/bin/env bash
#SBATCH --job-name=data-stage-230
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=executions/logs/data-stage-230-%j.out
#SBATCH --error=executions/logs/data-stage-230-%j.err

set -euo pipefail
: "${GSE105127_SOURCE_AUDIT_ROOT:?GSE105127_SOURCE_AUDIT_ROOT is required}"
: "${GSE105127_SOURCE_AUDIT_SHA256:?GSE105127_SOURCE_AUDIT_SHA256 is required}"
: "${GSE105127_V2_VALIDATION_ROOT:?GSE105127_V2_VALIDATION_ROOT is required}"
: "${GSE105127_V2_VALIDATION_SHA256:?GSE105127_V2_VALIDATION_SHA256 is required}"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src:${ROOT}"
"${PYTHON}" scripts/prepare_gse105127_gzip_preflight_revision_v2.py \
  --project-root "${ROOT}" \
  --source-campaign-root "${ROOT}/executions/gse105127-production-9675cd48993468d4" \
  --independent-audit-root "${GSE105127_SOURCE_AUDIT_ROOT}" \
  --independent-audit-sha256 "${GSE105127_SOURCE_AUDIT_SHA256}" \
  --validation-root "${GSE105127_V2_VALIDATION_ROOT}" \
  --validation-sha256 "${GSE105127_V2_VALIDATION_SHA256}" \
  --output-base "${ROOT}/executions"
