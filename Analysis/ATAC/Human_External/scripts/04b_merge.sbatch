#!/bin/bash
#SBATCH --job-name=snapmerge
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --output=Analysis/ATAC/Human_External/logs/merge_%j.out
#SBATCH --error=Analysis/ATAC/Human_External/logs/merge_%j.err
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
python Analysis/ATAC/Human_External/scripts/04b_merge_annotate.py
