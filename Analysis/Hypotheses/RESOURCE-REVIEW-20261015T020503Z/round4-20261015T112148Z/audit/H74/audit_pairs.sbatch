#!/bin/bash
#SBATCH --job-name=9742
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H74/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H74/logs/%j.err
set -euo pipefail
A=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H74
OUT=${A}/pairs-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}"
export TMPDIR=${A}/tmp XDG_CACHE_HOME=${A}/tmp PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
sha256sum GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz > ${OUT}/signal_pairs_sha256.txt
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python ${A}/audit_pairs.py > ${OUT}/audit_pairs.log 2>&1
sha256sum ${A}/audit_pairs.py ${A}/audit_pairs.sbatch > ${OUT}/script_sha256.tsv
