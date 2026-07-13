#!/bin/bash
#SBATCH --job-name=snapstep1
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --array=1-12
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=Analysis/ATAC/Human_External/logs/s1_%A_%a.out
#SBATCH --error=Analysis/ATAC/Human_External/logs/s1_%A_%a.err
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
python Analysis/ATAC/Human_External/scripts/04a_step1_one_donor.py
