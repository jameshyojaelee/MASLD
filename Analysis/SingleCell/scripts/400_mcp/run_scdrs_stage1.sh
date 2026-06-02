#!/bin/bash
#SBATCH --job-name=scdrs_s1_pp
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --output=logs/scdrs_stage1_%j.out
#SBATCH --error=logs/scdrs_stage1_%j.err

set -eo pipefail
mkdir -p logs

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

set -u
python Analysis/SingleCell/scripts/400_mcp/run_scdrs_stage1_preprocess.py
