#!/usr/bin/env bash
#SBATCH --job-name=data-stage-206
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=48:00:00
#SBATCH --output=executions/logs/data-stage-206-%j.out
#SBATCH --error=executions/logs/data-stage-206-%j.err

source slurm/gse105127_production_common.sh
if replay_complete "${CROSSWALK}"; then exit 0; fi
module purge
module load BEDTools/2.31.0-GCC-12.3.0 crossmap/0.7.0
begin_attempt crosswalk
STAGE=${ATTEMPT}/output.staging
"${PYTHON}" scripts/build_gse105127_cpg_crosswalk.py \
  --plan-root "${PLAN}" --reference-root "${REFERENCE}" \
  --collapsed-root "${RRBS_COLLAPSED}" --target-fasta "${TARGET_FASTA_PLAIN}" \
  --output "${STAGE}" --work "${ATTEMPT}/work"
publish_stage "${STAGE}" "${CROSSWALK}"
finish_attempt
