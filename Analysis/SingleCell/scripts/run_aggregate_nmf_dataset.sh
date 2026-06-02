#!/bin/bash
#SBATCH --job-name=nmf_agg_ds
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/nmf_agg_ds_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/nmf_agg_ds_%j.err

set -euo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT=$BASE
mkdir -p "$BASE/Analysis/SingleCell/scripts/logs"

cd "$BASE"
micromamba run -n spatial python "$BASE/Analysis/SingleCell/scripts/aggregate_nmf_celltype_stage_dataset.py"
