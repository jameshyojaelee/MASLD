#!/usr/bin/env bash
#SBATCH --job-name=data-stage-207
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-7
#SBATCH --output=executions/logs/data-stage-207-%A_%a.out
#SBATCH --error=executions/logs/data-stage-207-%A_%a.err

source slurm/gse105127_production_common.sh
module purge
module load RSEM/1.3.3-foss-2022a STAR/2.7.10b-GCC-11.3.0
begin_attempt rna-quant
"${PYTHON}" scripts/quantify_gse105127_rna_bundle.py \
  --plan-root "${PLAN}" --raw-root "${RNA_RAW}" \
  --reference-root "${RSEM_REFERENCE}" --bundle-id "${SLURM_ARRAY_TASK_ID}" \
  --output "${RNA_QUANT}" --attempt "${ATTEMPT}/members" \
  --threads "${SLURM_CPUS_PER_TASK}" --seed 105127
finish_attempt
