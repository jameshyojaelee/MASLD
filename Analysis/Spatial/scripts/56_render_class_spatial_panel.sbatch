#!/usr/bin/env bash
#SBATCH --job-name=figures
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/render_class_spatial_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/render_class_spatial_%j.out
set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
micromamba run -n rnaseq Rscript scripts/figures/fig4_hmsma_class_spatial.R
echo "EXIT=$?"
