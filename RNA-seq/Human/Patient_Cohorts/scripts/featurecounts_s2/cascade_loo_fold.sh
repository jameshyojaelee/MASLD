#!/bin/bash
#SBATCH --job-name=dream
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/loo_cv_s2/loo_%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/loo_cv_s2/loo_%x_%j.err

# Proper bash sbatch wrapper (NOT --wrap, which runs under /bin/sh where pipefail is illegal).
# HELD_OUT passed via --export. Partition + --mem set on the sbatch CLI per fold.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== LOO-CV -s 2 fold: HELD_OUT=${HELD_OUT} ($(date)) ==="
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/dream_loo_cv.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: dream_loo_cv.R HELD_OUT=${HELD_OUT}"; exit 1; fi
echo "Finished: $(date)"
