#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=48:00:00
#SBATCH --job-name=319_B2_reverse
#SBATCH --output=Analysis/SingleCell/logs/319_B2_reverse_%j.out
#SBATCH --error=Analysis/SingleCell/logs/319_B2_reverse_%j.err

set -eo pipefail
mkdir -p Analysis/SingleCell/logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript Analysis/SingleCell/scripts/319_nichenet_bulk_reverse_validation.R
