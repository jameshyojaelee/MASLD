#!/bin/bash
#SBATCH --job-name=6111
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H71-20261015T124239Z/logs/slurm-6111-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H71-20261015T124239Z/logs/slurm-6111-%j.err
set -euo pipefail
umask 027
H71=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H71-20261015T124239Z
STAMP=$(date -u +%Y%m%dT%H%M%SZ)
DA=${H71}/debug-${STAMP}/A_pre_outcome_real
DB=${H71}/debug-${STAMP}/B_null_lists
test ! -e "${H71}/debug-${STAMP}"; mkdir -p "${DA}/tmp" "${DB}/tmp"
export TMPDIR=${DA}/tmp XDG_CACHE_HOME=${DA}/tmp MPLCONFIGDIR=${DA}/tmp
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
ENVS=/gpfs/commons/home/jameslee/micromamba/envs
${ENVS}/spatial/bin/python ${H71}/scripts/01_xlsx_check.py > "${DA}/xlsx_check.out"
${ENVS}/rnaseq/bin/Rscript ${H71}/scripts/00_control_expression.R --out "${DA}"
/usr/bin/time -v ${ENVS}/rnaseq/bin/python ${H71}/scripts/02_h71_requalify.py --out "${DA}" --seed 20260923 --nperm 2000 --npcb 200 --nnc 100 --nnc2 200 --stop_before_outcome
cp "${DA}/xlsx_check.out" "${DA}/control_expression.tsv" "${DB}/"
/usr/bin/time -v ${ENVS}/rnaseq/bin/python ${H71}/scripts/02_h71_requalify.py --out "${DB}" --seed 20260923 --nperm 2000 --npcb 200 --nnc 100 --nnc2 200 --debug_null
