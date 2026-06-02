#!/bin/bash
#SBATCH --job-name=liver_11_12
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/scripts_11_12_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/scripts_11_12_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

echo "=== 11: Ortholog Mapping ==="
if [ $? -ne 0 ]; then echo "FAILED: 11"; exit 1; fi

echo ""
echo "=== 12: Library Intersection ==="
Rscript scripts/12_library_intersection.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 12"; exit 1; fi

echo ""
echo "=== Scripts 11-12 complete: $(date) ==="
