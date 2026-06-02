#!/bin/bash
#SBATCH --job-name=build_liver_markers
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/build_liver_markers_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/build_liver_markers_%j.err

set -e
mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs

export MAMBA_EXE="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MAMBA_EXE shell hook -s bash)"
micromamba activate spatial

python /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/build_liver_celltype_markers.py
