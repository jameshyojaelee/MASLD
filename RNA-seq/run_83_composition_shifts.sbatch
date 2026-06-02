#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=48:00:00
#SBATCH --job-name=83_comp_shift
#SBATCH --output=RNA-seq/logs/83_composition_%j.out
#SBATCH --error=RNA-seq/logs/83_composition_%j.err

set -eo pipefail
mkdir -p RNA-seq/logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript RNA-seq/83_composition_shifts.R
