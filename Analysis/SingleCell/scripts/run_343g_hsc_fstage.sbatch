#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --job-name=343g_hsc_fstage
#SBATCH --output=Analysis/SingleCell/logs/343g_hsc_fstage_%j.out
#SBATCH --error=Analysis/SingleCell/logs/343g_hsc_fstage_%j.err
set -eo pipefail
mkdir -p Analysis/SingleCell/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
micromamba run -n rapids_singlecell python -u Analysis/SingleCell/scripts/343g_hsc_activation_fstage.py
