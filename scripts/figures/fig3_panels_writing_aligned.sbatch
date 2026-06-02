#!/bin/bash
#SBATCH --job-name=fig3_panels
#SBATCH --partition=cpu
#SBATCH --mem=24G
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/fig3_panels_writing_aligned_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/fig3_panels_writing_aligned_%j.err

micromamba run -n rnaseq Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/fig3_panels_writing_aligned.R
