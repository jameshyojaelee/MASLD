#!/bin/bash
#SBATCH --job-name=aggreg_ct_stage
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=Analysis/SingleCell/logs/aggregate_celltype_stage_%j.out
#SBATCH --error=Analysis/SingleCell/logs/aggregate_celltype_stage_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

python Analysis/SingleCell/scripts/aggregate_celltype_stage_scores.py
