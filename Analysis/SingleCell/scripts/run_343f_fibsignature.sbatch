#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --job-name=343f_fibsig
#SBATCH --output=Analysis/SingleCell/logs/343f_fibsig_%j.out
#SBATCH --error=Analysis/SingleCell/logs/343f_fibsig_%j.err
set -eo pipefail
mkdir -p Analysis/SingleCell/logs
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba run -n rapids_singlecell python Analysis/SingleCell/scripts/343f_fibsignature_fstage.py
