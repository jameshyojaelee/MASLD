#!/bin/bash
#SBATCH --job-name=figure-context
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=figures/candidates/figure-context-%j.log
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -B scripts/manuscript/build_figure_claim_context_20261014.py
