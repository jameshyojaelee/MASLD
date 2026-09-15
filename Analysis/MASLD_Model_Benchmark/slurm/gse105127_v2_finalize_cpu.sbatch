#!/usr/bin/env bash
#SBATCH --job-name=data-stage-236
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=executions/logs/data-stage-236-%j.out
#SBATCH --error=executions/logs/data-stage-236-%j.err

source slurm/gse105127_gzip_v2_common.sh
if replay_complete "${TERMINAL}"; then exit 0; fi
begin_attempt finalize
if ! replay_complete "${ACTIVATED}"; then
  STAGE=${ATTEMPT}/activation.staging
  "${PYTHON}" scripts/finalize_gse105127_production_activation.py \
    --plan-root "${PLAN}" --reference-root "${REFERENCE}" \
    --collapsed-root "${RRBS_COLLAPSED}" --crosswalk-root "${CROSSWALK}" \
    --rna-root "${RNA_CONSOLIDATED}" --output "${STAGE}"
  publish_stage "${STAGE}" "${ACTIVATED}"
fi
STAGE=${ATTEMPT}/terminal.staging
"${PYTHON}" scripts/finalize_gse105127_gzip_preflight_revision_v2.py \
  --control-root "${CONTROL}" --activation-root "${ACTIVATED}" \
  --stage-links-root "${STAGE_LINKS}" \
  --independent-audit-root "${GSE105127_SOURCE_AUDIT_ROOT}" \
  --output "${STAGE}"
publish_stage "${STAGE}" "${TERMINAL}"
finish_attempt
