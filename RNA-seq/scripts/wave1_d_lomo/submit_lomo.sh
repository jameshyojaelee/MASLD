#!/bin/bash
#SBATCH --job-name=wave1d_lomo
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/wave1_d_lomo/lomo_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/wave1_d_lomo/lomo_%j.err

set -euo pipefail
source ~/.bashrc
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

Rscript RNA-seq/scripts/wave1_d_lomo/lomo_validation.R
