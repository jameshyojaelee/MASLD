#!/usr/bin/env bash
#SBATCH --job-name=data-stage-235
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=executions/logs/data-stage-235-%j.out
#SBATCH --error=executions/logs/data-stage-235-%j.err

source slurm/gse105127_gzip_v2_common.sh
if replay_complete "${RNA_CONSOLIDATED}"; then exit 0; fi
begin_attempt rna-consolidation
STAGE=${ATTEMPT}/output.staging
"${PYTHON}" scripts/consolidate_gse105127_rna_quantification.py \
  --plan-root "${PLAN}" --quant-root "${RNA_QUANT}" --output "${STAGE}"
publish_stage "${STAGE}" "${RNA_CONSOLIDATED}"
finish_attempt
