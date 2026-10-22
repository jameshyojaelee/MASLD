#!/bin/bash
#SBATCH --job-name=N3amb_model
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/N3-ambient-20261013T060951Z/%x-%j.out
# N3 ambient check, job B: aggregate to donors, unmodified N3 donor model on corrected counts, summaries.
# Usage: sbatch run_model.sbatch <STEP>   STEP = all | summary
set -euo pipefail
STEP="${1:-all}"
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RES=$REPO/Analysis/Hypotheses/N3-ambient-20261013T060951Z
N3=$REPO/Analysis/Hypotheses/N3-20261009T140326Z
SCR=$RES/scripts
N3SCR=$REPO/scripts/analysis/hypotheses_2026q4/N3_target_expression
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "host=$(hostname) job=$SLURM_JOB_ID step=$STEP start=$(date -u)"
(cd "$SCR" && sha256sum *.py *.R *.sh *.sbatch; sha256sum "$N3SCR/03_donor_model.R" "$RES/PRESPEC.md") > "$RES/script_sha256_job${SLURM_JOB_ID}.txt"
if [ "$STEP" = all ]; then
  "$RS" "$SCR/03_aggregate.R" "$RES" "$N3"
  "$RS" "$N3SCR/03_donor_model.R" "$RES" results_corrected 1000 "${SLURM_CPUS_PER_TASK:-8}"
fi
"$RS" "$SCR/04_ambient_summary.R" "$RES" "$N3"
echo "end=$(date -u)"
