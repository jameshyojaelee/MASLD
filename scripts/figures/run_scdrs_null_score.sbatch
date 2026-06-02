#!/bin/bash
#SBATCH --job-name=scdrs_null
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=700G
#SBATCH --cpus-per-task=4
#SBATCH --time=90:00:00
#SBATCH --array=0-199%4
#SBATCH --output=scripts/figures/logs/scdrs_null_%A_%a.out
#SBATCH --error=scripts/figures/logs/scdrs_null_%A_%a.err

# bigmem partition required: scDRS score_cell needs ~470GB peak (above cpu
# 210GB cap). interactive QoS used to bypass the 7TB nslab user-memory cap
# (other heavy jobs of this user already use most of the nslab quota);
# interactive caps us at 4 concurrent tasks, no memory ceiling.

set -eo pipefail
mkdir -p scripts/figures/logs

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial

set -u
export SCDRS_TRAIT_INDEX=$SLURM_ARRAY_TASK_ID
python scripts/figures/run_scdrs_null_score.py
