#!/bin/bash
#SBATCH --job-name=scdrs_s2_score
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=500G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --array=0-11
#SBATCH --output=logs/scdrs_stage2_%A_%a.out
#SBATCH --error=logs/scdrs_stage2_%A_%a.err

set -eo pipefail
mkdir -p logs

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

set -u
export SCDRS_TRAIT_INDEX=$SLURM_ARRAY_TASK_ID
python Analysis/SingleCell/scripts/400_mcp/run_scdrs_stage2_score.py
