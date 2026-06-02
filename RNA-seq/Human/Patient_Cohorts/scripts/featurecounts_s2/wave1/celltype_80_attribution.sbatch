#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=48:00:00
#SBATCH --job-name=TOAST
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/celltype_80_attribution_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/celltype_80_attribution_%j.err
# Wave1 [F]: 80 cell-type intrinsic attribution. Reads -s2 dream_results_ashr.csv + scRNA.
# Writes celltype_primary_attribution.csv -> GATE for 85/87/88. NOT an atlas feeder.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== 80 celltype intrinsic attribution | $(date) ==="
Rscript RNA-seq/80_celltype_intrinsic_attribution.R
if [ $? -ne 0 ]; then echo "FAILED: 80_celltype_intrinsic_attribution"; exit 1; fi
echo "=== done $(date) ==="
