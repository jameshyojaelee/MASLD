#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=2 --mem=10G --time=48:00:00
#SBATCH --job-name=30_I1 --output=Analysis/Spatial/logs/30_I1_%j.out --error=Analysis/Spatial/logs/30_I1_%j.err
set -eo pipefail; mkdir -p Analysis/Spatial/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript Analysis/Spatial/scripts/30_spatial_ccc_consensus.R
