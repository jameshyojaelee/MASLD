#!/usr/bin/env bash
#SBATCH --job-name=data-stage-204
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --array=0-7
#SBATCH --output=executions/logs/data-stage-204-%A_%a.out
#SBATCH --error=executions/logs/data-stage-204-%A_%a.err

source slurm/gse105127_production_common.sh
begin_attempt rna-download
"${PYTHON}" scripts/download_gse105127_rna_bundle.py \
  --plan-root "${PLAN}" --bundle-id "${SLURM_ARRAY_TASK_ID}" \
  --output "${RNA_RAW}" --attempt "${ATTEMPT}/members"
finish_attempt
