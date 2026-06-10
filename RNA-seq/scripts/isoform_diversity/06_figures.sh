#!/bin/bash
#SBATCH --job-name=isofigures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=logs/figures_%x_%j.out
#SBATCH --error=logs/figures_%x_%j.err
# sbatch --export=ALL,SPECIES=human 06_figures.sh   (runs in rnaseq env for publication_theme)
set -euo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ/RNA-seq/scripts/isoform_diversity"
micromamba run -n rnaseq Rscript 06_figures.R "${SPECIES:-human}"
