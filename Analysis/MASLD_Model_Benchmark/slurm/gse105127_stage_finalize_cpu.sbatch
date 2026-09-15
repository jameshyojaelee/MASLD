#!/usr/bin/env bash
#SBATCH --job-name=data-stage-209
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=executions/logs/data-stage-209-%j.out
#SBATCH --error=executions/logs/data-stage-209-%j.err

source slurm/gse105127_production_common.sh
if replay_complete "${ACTIVATED}"; then exit 0; fi
begin_attempt finalize
STAGE=${ATTEMPT}/output.staging
"${PYTHON}" scripts/finalize_gse105127_production_activation.py \
  --plan-root "${PLAN}" --reference-root "${REFERENCE}" \
  --collapsed-root "${RRBS_COLLAPSED}" --crosswalk-root "${CROSSWALK}" \
  --rna-root "${RNA_CONSOLIDATED}" --output "${STAGE}"
publish_stage "${STAGE}" "${ACTIVATED}"
finish_attempt
