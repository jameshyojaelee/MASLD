#!/bin/bash
#SBATCH --job-name=ggplot
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/proteo_cascade_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/proteo_cascade_%j.err
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript scripts/figures/proteo_transcript_cascade.R
echo "PROTEO_CASCADE_DONE $(date)"
