#!/bin/bash
#SBATCH --job-name=348_fstage_pb
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=80G
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/scripts/logs/348_pseudobulk_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/logs/348_pseudobulk_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

source /gpfs/commons/home/jameslee/micromamba/etc/profile.d/micromamba.sh
micromamba activate spatial

python Analysis/SingleCell/scripts/348_celltype_fstage_pseudobulk.py
