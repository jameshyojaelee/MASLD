#!/bin/bash
#SBATCH --job-name=dream
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=%x_lvqw04_%j.out
#SBATCH --error=%x_lvqw04_%j.out
# ---------------------------------------------------------------------------
# LVQW staging run — Module 04 (limma-voom-quality-weighted, dataset fixed).
# Writes to lvqw_staging/ ONLY (SEXV3_OUT override). Canonical sex_v3 untouched.
# ---------------------------------------------------------------------------
set -euo pipefail
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
export SEXV3_OUT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3/lvqw_staging
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

cd "$SCRIPTS"
source /gpfs/commons/home/jameslee/micromamba/etc/profile.d/micromamba.sh
micromamba activate rnaseq

echo "=== LVQW staging Module 04 ==="
echo "host: $(hostname)  cpus: ${SLURM_CPUS_PER_TASK:-NA}  SEXV3_OUT: $SEXV3_OUT"
echo "started: $(date)"
Rscript 04_dream_M2_extended.R
echo "finished 04: $(date)"
