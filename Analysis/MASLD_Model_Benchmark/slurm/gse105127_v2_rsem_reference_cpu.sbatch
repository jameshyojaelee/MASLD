#!/usr/bin/env bash
#SBATCH --job-name=data-stage-231
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=executions/logs/data-stage-231-%j.out
#SBATCH --error=executions/logs/data-stage-231-%j.err

source slurm/gse105127_gzip_v2_common.sh
LINK=${STAGE_LINKS}/rsem_reference
if replay_complete "${LINK}"; then exit 0; fi
module purge
module load RSEM/1.3.3-foss-2022a STAR/2.7.10b-GCC-11.3.0
begin_attempt rsem-reference
run_preflight rsem_reference rsem_reference \
  --target-fasta "${TARGET_FASTA}" --target-gtf "${TARGET_GTF}"
if ! replay_complete "${RSEM_REFERENCE}"; then
  STAGE=${ATTEMPT}/output.staging
  "${PYTHON}" scripts/build_gse105127_rsem_reference.py \
    --plan-root "${PLAN}" --target-fasta "${TARGET_FASTA}" \
    --target-gtf "${TARGET_GTF}" --output "${STAGE}" \
    --threads "${SLURM_CPUS_PER_TASK}"
  publish_stage "${STAGE}" "${RSEM_REFERENCE}"
fi
freeze_stage_link rsem_reference "${RSEM_REFERENCE}" "${LINK}"
finish_attempt
