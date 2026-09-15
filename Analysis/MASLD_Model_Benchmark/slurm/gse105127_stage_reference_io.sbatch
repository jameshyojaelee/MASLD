#!/usr/bin/env bash
#SBATCH --job-name=data-stage-202
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=executions/logs/data-stage-202-%j.out
#SBATCH --error=executions/logs/data-stage-202-%j.err

source slurm/gse105127_production_common.sh
if replay_complete "${REFERENCE}"; then exit 0; fi
module purge
module load SAMtools/1.18-GCC-12.3.0
begin_attempt reference
STAGE=${ATTEMPT}/output.staging
"${PYTHON}" scripts/freeze_gse105127_reference_bundle.py \
  --plan-root "${PLAN}" --target-fasta "${TARGET_FASTA}" \
  --target-gtf "${TARGET_GTF}" --output "${STAGE}"
publish_stage "${STAGE}" "${REFERENCE}"
finish_attempt
