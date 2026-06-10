#!/bin/bash
#SBATCH --job-name=dream
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=%x_lvqw0506_%j.out
#SBATCH --error=%x_lvqw0506_%j.out
# ---------------------------------------------------------------------------
# LVQW staging run — Module 05 (mashr) then 06 (classify), chained after 04.
# Reads + writes lvqw_staging/ ONLY (SEXV3_OUT override). V2_F/CHRX_F canonical.
# ---------------------------------------------------------------------------
set -euo pipefail
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
export SEXV3_OUT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3/lvqw_staging
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

cd "$SCRIPTS"
source /gpfs/commons/home/jameslee/micromamba/etc/profile.d/micromamba.sh
micromamba activate rnaseq

echo "=== LVQW staging Module 05 (mashr) ==="
echo "host: $(hostname)  SEXV3_OUT: $SEXV3_OUT"
echo "started 05: $(date)"
Rscript 05_mashr_fit.R
echo "finished 05: $(date)"

echo "=== LVQW staging Module 06 (classify) ==="
echo "started 06: $(date)"
Rscript 06_classify_v3.R
echo "finished 06: $(date)"
