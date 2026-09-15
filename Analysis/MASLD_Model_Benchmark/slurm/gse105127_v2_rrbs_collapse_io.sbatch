#!/usr/bin/env bash
#SBATCH --job-name=data-stage-232
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=20:00:00
#SBATCH --array=0-7
#SBATCH --output=executions/logs/data-stage-232-%A_%a.out
#SBATCH --error=executions/logs/data-stage-232-%A_%a.err

source slurm/gse105127_gzip_v2_common.sh
LINK=${STAGE_LINKS}/rrbs_collapse_bundle_$(printf '%02d' "${SLURM_ARRAY_TASK_ID}")
if replay_complete "${LINK}"; then exit 0; fi
begin_attempt rrbs-collapse
run_preflight rrbs_collapse rrbs_collapse \
  --plan-root "${PLAN}" --bundle-id "${SLURM_ARRAY_TASK_ID}"
"${PYTHON}" scripts/collapse_gse105127_rrbs_bundle.py \
  --plan-root "${PLAN}" --bundle-id "${SLURM_ARRAY_TASK_ID}" \
  --output "${RRBS_COLLAPSED}" --attempt "${ATTEMPT}/members"
freeze_stage_link rrbs_collapse "${RRBS_COLLAPSED}" "${LINK}" \
  --bundle-id "${SLURM_ARRAY_TASK_ID}"
finish_attempt
