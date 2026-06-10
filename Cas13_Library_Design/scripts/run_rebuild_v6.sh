#!/bin/bash
#SBATCH --job-name=cas13lib
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=Cas13_Library_Design/logs/rebuild_v6_%j.out
#SBATCH --error=Cas13_Library_Design/logs/rebuild_v6_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

Rscript Cas13_Library_Design/scripts/rebuild_cas13_library.R
