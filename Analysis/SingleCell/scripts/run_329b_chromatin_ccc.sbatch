#!/bin/bash
#SBATCH --job-name=chromatin_ccc_v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/329b_chromatin_ccc_%j.out
#SBATCH --error=logs/329b_chromatin_ccc_%j.err

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell
mkdir -p logs results_gpu_v2/chromatin_ccc
Rscript scripts/329b_chromatin_ccc_celltype_resolved.R
