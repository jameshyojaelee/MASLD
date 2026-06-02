#!/bin/bash
#SBATCH --job-name=GSE246088_DE
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets/GSE246088/scripts/logs/de_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets/GSE246088/scripts/logs/de_%j.err

echo "=== GSE246088 DE Analysis ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Date: $(date)"

# Activate rnaseq environment (set -u deferred past conda activate)
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -euo pipefail

Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets/GSE246088/scripts/run_de_analysis.R

echo "=== DONE ==="
echo "Date: $(date)"
