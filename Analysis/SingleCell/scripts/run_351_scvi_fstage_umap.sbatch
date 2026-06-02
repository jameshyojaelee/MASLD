#!/bin/bash
#SBATCH --job-name=351_scvi_fstage_umap
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/351_scvi_fstage_umap_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/351_scvi_fstage_umap_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p Analysis/SingleCell/scripts/logs

eval "$(micromamba shell hook --shell=bash)"
micromamba activate rapids_singlecell

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

python -u Analysis/SingleCell/scripts/351_scvi_fstage_umap_figure.py
