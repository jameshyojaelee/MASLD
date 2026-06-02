#!/bin/bash
#SBATCH --job-name=dream
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/s2_bariatric_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/s2_bariatric_%j.err

set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

# STAR-native bariatric sensitivity: re-run dream excluding GSE162694 (Bril bariatric cohort).
# Script 05 reads canonical -s 2 merged_dge.rds; EXCL_EXTRA drops the cohort and
# writes dream_results_exclGSE162694.csv (compare to canonical for Jaccard/rho).
export EXCL_EXTRA="GSE162694"

echo "=== Bariatric sensitivity (STAR -s 2; exclude GSE162694) ==="
echo "Started: $(date)"
Rscript analysis/integration/scripts/05_dream_mega_analysis.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: bariatric dream"; exit 1; fi
echo "Finished: $(date)"
