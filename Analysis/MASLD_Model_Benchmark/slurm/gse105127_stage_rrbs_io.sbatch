#!/usr/bin/env bash
#SBATCH --job-name=data-stage-203
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=20:00:00
#SBATCH --array=0-7
#SBATCH --output=executions/logs/data-stage-203-%A_%a.out
#SBATCH --error=executions/logs/data-stage-203-%A_%a.err

source slurm/gse105127_production_common.sh
begin_attempt rrbs
"${PYTHON}" scripts/collapse_gse105127_rrbs_bundle.py \
  --plan-root "${PLAN}" --bundle-id "${SLURM_ARRAY_TASK_ID}" \
  --output "${RRBS_COLLAPSED}" --attempt "${ATTEMPT}/members"
finish_attempt
