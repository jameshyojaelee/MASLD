#!/bin/bash -l
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=B3_xtab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/outputs/team_B/B3_logs/xtab_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/outputs/team_B/B3_logs/xtab_%j.err

# B3 sex x ancestry x COLOC cross-tab.
# Consumes existing pooled 28-GWAS COLOC + (when available) Pan-UKBB
# and FinnGen sex-stratified COLOC. See 216b_sex_ancestry_coloc_xtab.R.

set -eo pipefail

BASE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE_DIR}/RNA-seq"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "=== B3 Sex x Ancestry x COLOC cross-tab ==="
echo "Job: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo

Rscript 216b_sex_ancestry_coloc_xtab.R

echo
echo "End: $(date)"
