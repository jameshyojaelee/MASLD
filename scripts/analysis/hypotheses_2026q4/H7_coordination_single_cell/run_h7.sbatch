#!/bin/bash
#SBATCH --job-name=H7count
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
# Usage: sbatch run_h7.sbatch <OUT_DIR>
set -euo pipefail
OUT="${1:?out dir required}"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript \
  /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H7_coordination_single_cell/count_f0f1_donors.R "$OUT"
