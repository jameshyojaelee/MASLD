#!/bin/bash
#SBATCH --job-name=6113
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H71-20261015T124239Z/logs/slurm-6113-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H71-20261015T124239Z/logs/slurm-6113-%j.err
set -euo pipefail
umask 027
H71=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H71-20261015T124239Z
DB=${H71}/debug-$(date -u +%Y%m%dT%H%M%SZ)/B_null_lists
test ! -e "${DB}"; mkdir -p "${DB}/tmp"
export TMPDIR=${DB}/tmp XDG_CACHE_HOME=${DB}/tmp MPLCONFIGDIR=${DB}/tmp
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
ENVS=/gpfs/commons/home/jameslee/micromamba/envs
cp "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H71-20261015T124239Z/debug-20261015T125337Z/A_pre_outcome_real/xlsx_check.out" "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H71-20261015T124239Z/debug-20261015T125337Z/A_pre_outcome_real/control_expression.tsv" "${DB}/"
${ENVS}/rnaseq/bin/python ${H71}/scripts/02_h71_requalify.py --out "${DB}" --seed 20260923 --nperm 2000 --npcb 200 --nnc 100 --nnc2 200 --debug_null
