#!/usr/bin/env bash
#SBATCH --job-name=data-stage-214
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=executions/logs/data-stage-214-%j.out
#SBATCH --error=executions/logs/data-stage-214-%j.err

source slurm/gse105127_production_common.sh
if replay_complete "${REFERENCE}"; then exit 0; fi
module purge
module load SAMtools/1.18-GCC-12.3.0
EXACT_SCRIPT=scripts/freeze_gse105127_reference_bundle_legacy_exact.py
EXACT_SCRIPT_SHA256=e3610b025be3c4c9238fb111c4225bdde47f69eb22cc9ebe83309eff3b8b462e
printf '%s  %s\n' "${EXACT_SCRIPT_SHA256}" "${EXACT_SCRIPT}" | sha256sum --check --strict
begin_attempt reference-exact
STAGE=${ATTEMPT}/output.staging
"${PYTHON}" "${EXACT_SCRIPT}" \
  --plan-root "${PLAN}" \
  --target-fasta "${TARGET_FASTA}" \
  --target-gtf "${TARGET_GTF}" \
  --legacy-diagnostic-root "${ROOT}/executions/data-check-211-21084942" \
  --alphabet-diagnostic-root "${ROOT}/executions/data-check-216-21086818" \
  --source-cache-root "${ROOT}/executions/gse105127-reference-source-cache-21087810" \
  --implementation-sha256 "${EXACT_SCRIPT_SHA256}" \
  --output "${STAGE}"
publish_stage "${STAGE}" "${REFERENCE}"
finish_attempt
