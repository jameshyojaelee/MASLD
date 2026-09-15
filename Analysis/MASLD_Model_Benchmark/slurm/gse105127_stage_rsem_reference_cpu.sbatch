#!/usr/bin/env bash
#SBATCH --job-name=data-stage-205
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=executions/logs/data-stage-205-%j.out
#SBATCH --error=executions/logs/data-stage-205-%j.err

source slurm/gse105127_production_common.sh
if replay_complete "${RSEM_REFERENCE}"; then exit 0; fi
module purge
module load RSEM/1.3.3-foss-2022a STAR/2.7.10b-GCC-11.3.0
begin_attempt rsem-reference
STAGE=${ATTEMPT}/output.staging
"${PYTHON}" scripts/build_gse105127_rsem_reference.py \
  --plan-root "${PLAN}" --target-fasta "${TARGET_FASTA}" \
  --target-gtf "${TARGET_GTF}" --output "${STAGE}" \
  --threads "${SLURM_CPUS_PER_TASK}"
publish_stage "${STAGE}" "${RSEM_REFERENCE}"
finish_attempt
