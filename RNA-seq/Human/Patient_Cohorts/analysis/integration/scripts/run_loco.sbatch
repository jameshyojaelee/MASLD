#!/bin/bash
#SBATCH --job-name=LOCO
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --output=logs/LOCO_%x_%j.out
#SBATCH --error=logs/LOCO_%x_%j.err

# LOO-CV dream refit: holds out cohort specified by HELD_OUT env var.
# Submit once per cohort:
#   for c in GSE126848 GSE130970 GSE135251 GSE162694 GSE213621; do
#       sbatch --export=ALL,HELD_OUT=$c run_loco.sbatch
#   done

set -eo pipefail
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

if [ -z "${HELD_OUT:-}" ]; then
    echo "ERROR: HELD_OUT env var must be set (e.g., GSE213621)" >&2
    exit 1
fi

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
echo "[$(date)] LOO-CV: HELD_OUT=${HELD_OUT} on $(hostname)"
Rscript dream_loo_cv.R
echo "[$(date)] Done LOO-CV: ${HELD_OUT}"
