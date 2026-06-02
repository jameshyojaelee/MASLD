#!/bin/bash
#SBATCH --job-name=348b_dream
#SBATCH --partition=cpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/scripts/logs/348b_dream_%A_%a.out
#SBATCH --error=Analysis/SingleCell/scripts/logs/348b_dream_%A_%a.err
#SBATCH --array=0-1

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

source /gpfs/commons/home/jameslee/micromamba/etc/profile.d/micromamba.sh
micromamba activate rnaseq

TAGS=(augmented documented)
TAG=${TAGS[$SLURM_ARRAY_TASK_ID]}
echo "[run] COHORT_TAG=$TAG"

COHORT_TAG=$TAG Rscript Analysis/SingleCell/scripts/348b_celltype_fstage_dream.R
