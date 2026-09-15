#!/usr/bin/env bash
#SBATCH --job-name=data-stage-233
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-7
#SBATCH --output=executions/logs/data-stage-233-%A_%a.out
#SBATCH --error=executions/logs/data-stage-233-%A_%a.err

source slurm/gse105127_gzip_v2_common.sh
LINK=${STAGE_LINKS}/rna_quantification_bundle_$(printf '%02d' "${SLURM_ARRAY_TASK_ID}")
if replay_complete "${LINK}"; then exit 0; fi
module purge
module load RSEM/1.3.3-foss-2022a STAR/2.7.10b-GCC-11.3.0
begin_attempt rna-quantification
run_preflight rna_quantification rna_quantification \
  --plan-root "${PLAN}" --raw-root "${RNA_RAW}" \
  --bundle-id "${SLURM_ARRAY_TASK_ID}"
"${PYTHON}" scripts/quantify_gse105127_rna_bundle_v2.py \
  --plan-root "${PLAN}" --raw-root "${RNA_RAW}" \
  --reference-root "${RSEM_REFERENCE}" --bundle-id "${SLURM_ARRAY_TASK_ID}" \
  --output "${RNA_QUANT}" --attempt "${ATTEMPT}/members" \
  --threads "${SLURM_CPUS_PER_TASK}" --seed 105127
freeze_stage_link rna_quantification "${RNA_QUANT}" "${LINK}" \
  --bundle-id "${SLURM_ARRAY_TASK_ID}"
finish_attempt
