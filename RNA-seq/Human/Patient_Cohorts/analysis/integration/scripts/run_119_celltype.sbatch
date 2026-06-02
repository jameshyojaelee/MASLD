#!/bin/bash
#SBATCH --job-name=119_celltype_prog
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/119_celltype_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/119_celltype_%j.err

set -euo pipefail

eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq

echo "Running Script 119 on $(hostname) at $(date)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/119_celltype_resolved_progression.R

echo "Finished at $(date)"
