#!/usr/bin/env bash
#SBATCH --job-name=data-stage-201
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=executions/logs/data-stage-201-%j.out
#SBATCH --error=executions/logs/data-stage-201-%j.err

source slurm/gse105127_production_common.sh
if replay_complete "${PLAN}"; then exit 0; fi
begin_attempt plan
STAGE=${ATTEMPT}/output.staging
"${PYTHON}" scripts/plan_gse105127_production_dag.py \
  --activation-root "${ACTIVATION}" --rrbs-root "${RRBS_SOURCE}" \
  --output "${STAGE}" --bundles 8
publish_stage "${STAGE}" "${PLAN}"
finish_attempt
