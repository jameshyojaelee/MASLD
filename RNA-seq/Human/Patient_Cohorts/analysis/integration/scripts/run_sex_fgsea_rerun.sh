#!/bin/bash
#SBATCH --job-name=sex_fgsea_rerun
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/sex_fgsea_rerun_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/sex_fgsea_rerun_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00

# Rerun script 26 to fix sex-specific fgsea (switch from logFC to t-stat ranking)
set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

echo "=== Sex fgsea rerun: $(date) ==="
echo "Fix: t-statistic ranking + nPermSimple=10000 (matching stage GSEA scripts)"
echo ""

Rscript "$SCRIPTS/26_sex_stratified_analysis.R" 2>&1

echo ""
echo "=== Completed: $(date) ==="
